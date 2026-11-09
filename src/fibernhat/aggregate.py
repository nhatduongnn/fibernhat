"""Aggregate (piled-up) fiber data around a set of sites, straight from BAM.

Pure numpy replacement of the pyranges path
(Fiber_seq_plots.process_aggregate_fiber_plot + aggregate_fiberplot_open_sites
WindowedDecodeMixin): every read is decoded once, each (read, site) overlap
becomes one row holding the oriented positions (distance to site, site
5'->3' to the right) of its methylated and called A's.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .io.bam import ReadFilter, Thresholds, open_bam, read_label, read_to_fiber


@dataclass
class SiteFibers:
    half: int
    labels: List[str]
    rows: pd.DataFrame            # site_idx, label, read_id, strand, start, end, n_meth_win, n_called_win
    meth_ptr: np.ndarray
    meth_rel: np.ndarray          # int32 oriented positions (-half..half)
    called_ptr: np.ndarray
    called_rel: np.ndarray
    mod: Dict[str, np.ndarray]    # per label, length 2*half+1
    valid: Dict[str, np.ndarray]
    n_sites_used: int

    def meth(self, i: int) -> np.ndarray:
        return self.meth_rel[self.meth_ptr[i]:self.meth_ptr[i + 1]]

    def profile(self, label: str) -> np.ndarray:
        with np.errstate(all="ignore"):
            return self.mod[label] / self.valid[label]

    @property
    def x(self) -> np.ndarray:
        return np.arange(-self.half, self.half + 1)

    def save(self, path: str):
        np.savez_compressed(path, half=self.half, labels=np.array(self.labels),
                            row_cols=np.array(list(self.rows.columns)),
                            **{f"row__{c}": (np.array(self.rows[c].astype(str).tolist(), dtype="U")
                                             if not pd.api.types.is_numeric_dtype(self.rows[c]) else self.rows[c].values)
                               for c in self.rows.columns}, meth_ptr=self.meth_ptr,
                            meth_rel=self.meth_rel, called_ptr=self.called_ptr,
                            called_rel=self.called_rel, n_sites_used=self.n_sites_used,
                            mod=np.array([self.mod[l] for l in self.labels]),
                            valid=np.array([self.valid[l] for l in self.labels]))

    @classmethod
    def load(cls, path: str) -> "SiteFibers":
        z = np.load(path, allow_pickle=False)
        labels = [str(l) for l in z["labels"]]
        return cls(int(z["half"]), labels, pd.DataFrame({str(c): z[f"row__{c}"] for c in z["row_cols"]}), z["meth_ptr"], z["meth_rel"],
                   z["called_ptr"], z["called_rel"],
                   dict(zip(labels, z["mod"])), dict(zip(labels, z["valid"])), int(z["n_sites_used"]))


def collect(bam_paths: Sequence[str], sites: pd.DataFrame, half: int, thresholds: Thresholds,
            filt: Optional[ReadFilter] = None, label_by: str = "auto", require_span: bool = True,
            max_per_site: Optional[int] = None, sort_win: int = 100, seed: int = 0,
            chrom_lens: Optional[Dict[str, int]] = None) -> SiteFibers:
    filt = filt if filt is not None else ReadFilter()
    W = 2 * half + 1
    sites = sites.reset_index(drop=True)
    rows, meth_parts, called_parts = [], [], []
    mod: Dict[str, np.ndarray] = {}
    valid: Dict[str, np.ndarray] = {}
    used_sites = set()
    for p in bam_paths:
        with open_bam(p) as bam:
            refs = set(bam.references)
            lens = dict(zip(bam.references, bam.lengths))
            for chrom, g in sites.groupby("chrom", sort=False):
                if chrom not in refs:
                    continue
                g = g[(g.center - half >= 0) & (g.center + half < lens[chrom])].sort_values("center")
                if g.empty:
                    continue
                centers = g.center.values
                signs = g.sign.values
                sidx = g.index.values
                lo, hi = int(centers.min() - half), int(centers.max() + half + 1)
                for read in bam.fetch(chrom, max(0, lo), hi):
                    if filt.reason(read) is not None:
                        continue
                    rs, re_ = read.reference_start, read.reference_end
                    if require_span:
                        a = np.searchsorted(centers, rs + half, "left")
                        b = np.searchsorted(centers, re_ - half - 1, "right")
                    else:
                        a = np.searchsorted(centers, rs - half, "left")
                        b = np.searchsorted(centers, re_ + half, "right")
                    if b <= a:
                        continue
                    label = read_label(read, p, label_by)
                    f = read_to_fiber(read, thresholds.get(label), label)
                    called = f.call >= 0
                    cpos, cmeth = f.pos[called], f.call[called] == 1
                    if label not in mod:
                        mod[label] = np.zeros(W, np.int64)
                        valid[label] = np.zeros(W, np.int64)
                    for j in range(a, b):
                        c, s = int(centers[j]), int(signs[j])
                        lo_i = np.searchsorted(cpos, c - half, "left")
                        hi_i = np.searchsorted(cpos, c + half, "right")
                        rel = (cpos[lo_i:hi_i] - c) * s
                        m = cmeth[lo_i:hi_i]
                        if rel.size == 0:
                            continue
                        np.add.at(valid[label], rel + half, 1)
                        np.add.at(mod[label], rel[m] + half, 1)
                        mr = np.sort(rel[m]).astype(np.int32)
                        meth_parts.append(mr)
                        called_parts.append(np.sort(rel).astype(np.int32))
                        inwin = np.abs(rel) <= sort_win
                        rows.append((int(sidx[j]), label, f.read_id, f.strand, rs, re_,
                                     int((m & inwin).sum()), int(inwin.sum())))
                        used_sites.add(int(sidx[j]))
    rows_df = pd.DataFrame(rows, columns=["site_idx", "label", "read_id", "strand", "start", "end",
                                          "n_meth_win", "n_called_win"])
    if max_per_site is not None and len(rows_df):
        rng = np.random.default_rng(seed)
        rows_df["_r"] = rng.random(len(rows_df))
        keep = rows_df.sort_values("_r").groupby("site_idx").head(max_per_site).index.sort_values()
        rows_df = rows_df.loc[keep].drop(columns="_r")
        meth_parts = [meth_parts[i] for i in keep]
        called_parts = [called_parts[i] for i in keep]
        rows_df = rows_df.reset_index(drop=True)

    def csr(parts):
        ptr = np.zeros(len(parts) + 1, np.int64)
        ptr[1:] = np.cumsum([len(x) for x in parts])
        return ptr, (np.concatenate(parts) if parts else np.zeros(0, np.int32))

    mp, mr = csr(meth_parts)
    cp, cr = csr(called_parts)
    if rows_df.size:
        rows_df = rows_df.merge(sites[["chrom", "center", "sign", "name"]], left_on="site_idx",
                                right_index=True, how="left")
    return SiteFibers(half, sorted(mod), rows_df, mp, mr, cp, cr, mod, valid, len(used_sites))


def order_rows(sf: SiteFibers, idx: np.ndarray, how: str = "shuffle", seed: int = 0) -> np.ndarray:
    """Order a subset of row indices. how: shuffle | meth_in_win | meth_frac | site | none."""
    r = sf.rows.iloc[idx]
    if how == "shuffle":
        return np.random.default_rng(seed).permutation(idx)
    if how == "meth_in_win":
        frac = r.n_meth_win / r.n_called_win.clip(lower=1)
        return idx[np.argsort(-frac.values, kind="stable")]
    if how == "meth_frac":
        fr = np.array([len(sf.meth(i)) / max(1, sf.called_ptr[i + 1] - sf.called_ptr[i]) for i in idx])
        return idx[np.argsort(-fr, kind="stable")]
    if how == "site":
        return idx[np.lexsort((r.start.values, r.site_idx.values))]
    return idx


def profile_table(sf: SiteFibers, site_name: str) -> pd.DataFrame:
    out = []
    for lab in sf.labels:
        out.append(pd.DataFrame({"site_set": site_name, "label": lab, "pos": sf.x,
                                 "n_mod": sf.mod[lab], "n_valid": sf.valid[lab],
                                 "m6a_frac": sf.profile(lab)}))
    return pd.concat(out, ignore_index=True)
