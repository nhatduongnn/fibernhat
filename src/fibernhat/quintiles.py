"""Per-read m6A-rate quintiles and their metaprofiles at dyads / TSS / TFs.

Port of investigate_un_crosslinked_runs step 6 (per_read_heterogeneity.py +
footprint_metaprofiles.py), BAM-only:

pass 1  per read (filter passed, nuclear, read length >= min_len):
        rate = methylated / called target bases over the whole stored read
        (soft clips included, as in the original); quintiles per label with
        pd.qcut(rate, 5) -> Q1 = least methylated ... Q5 = most.
pass 2  per (label, quintile): per-position mod/valid pileup (both strands
        summed, engine pysam or modkit), summed over oriented site windows;
        plus the genome-wide (nuclear) m6A rate of that quintile.

Original settings: ML>=128 binary call, MAPQ>=10, no qs/rDNA filter.
Reproduce with --threshold-ml 128 --min-mapq 10 --min-qs none --keep-rdna.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .io import refs
from .io.bam import ReadFilter, Thresholds, iter_reads, open_bam, read_length, read_space_counts, read_to_fiber
from . import pileup as P

N_Q = 5


def per_read_rates(bam_paths: Sequence[str], thresholds: Thresholds, filt: ReadFilter,
                   label_by: str = "auto", min_len: int = 500) -> pd.DataFrame:
    rows = []
    for read, label, _ in iter_reads(bam_paths, None, filt, label_by):
        if read.reference_name not in refs.NUCLEAR:
            continue
        L = read_length(read)
        if L < min_len:
            continue
        m, c, n = read_space_counts(read, thresholds.get(label))
        rows.append((read.query_name, label, read.reference_name, read.reference_start,
                     read.reference_end, L, n, m, c))
    df = pd.DataFrame(rows, columns=["read_id", "label", "chrom", "start", "end", "read_len",
                                     "n_target", "n_meth", "n_canon"])
    df["rate"] = df.n_meth / (df.n_meth + df.n_canon).clip(lower=1)
    df["quintile"] = 0
    for lab, g in df.groupby("label"):
        if len(g) >= N_Q:
            df.loc[g.index, "quintile"] = pd.qcut(g.rate.rank(method="first"), N_Q, labels=False) + 1
    return df


def window_sums(mod_valid: Dict[str, np.ndarray], sites: pd.DataFrame, half: int):
    """mod_valid[chrom] = int[2, L] (mod, valid). Returns (mod[W], valid[W], n_used)."""
    W = 2 * half + 1
    tm, tv, n = np.zeros(W), np.zeros(W), 0
    for chrom, g in sites.groupby("chrom"):
        a = mod_valid.get(chrom)
        if a is None:
            continue
        L = a.shape[1]
        for c, s in zip(g.center.values, g.sign.values):
            lo, hi = c - half, c + half + 1
            if lo < 0 or hi > L:
                continue
            tm += a[0, lo:hi][::s]
            tv += a[1, lo:hi][::s]
            n += 1
    return tm, tv, n


def _strand_sum(arrays: P.Arrays) -> Dict[str, np.ndarray]:
    return {c: np.stack([a[0] + a[2], a[1] + a[3]]) for c, a in arrays.items()}


def run(bam_paths: Sequence[str], thresholds: Thresholds, filt: ReadFilter,
        site_specs: Sequence[str] = ("dyad", "tss", "abf1", "reb1"),
        halves: Optional[Dict[str, int]] = None, label_by: str = "auto", min_len: int = 500,
        engine: str = "pysam", modkit: Optional[str] = None, threads: int = 4,
        progress: bool = True):
    """Returns (reads DataFrame, profiles DataFrame, summary DataFrame)."""
    import sys
    reads = per_read_rates(bam_paths, thresholds, filt, label_by, min_len)
    if progress:
        print(f"[quintiles] pass 1: {len(reads):,} reads", file=sys.stderr)
    sites = {s: refs.load_sites(s) for s in site_specs}
    halves = {s: (halves or {}).get(s, refs.DEFAULT_HALF.get(s, 1000)) for s in site_specs}
    prof_rows, summ = [], []
    keys = [(lab, q) for lab in sorted(reads.label.unique()) for q in range(1, N_Q + 1)
            if ((reads.label == lab) & (reads.quintile == q)).any()]
    W = {s: 2 * halves[s] + 1 for s in site_specs}
    acc = {k: {"gm": 0, "gv": 0, **{s: [np.zeros(W[s]), np.zeros(W[s]), 0] for s in site_specs}}
           for k in keys}

    def consume(k, chrom, mv2):
        a = acc[k]
        if chrom in refs.NUCLEAR:
            a["gm"] += int(mv2[0].sum()); a["gv"] += int(mv2[1].sum())
        for s in site_specs:
            tm, tv, n = window_sums({chrom: mv2}, sites[s], halves[s])
            a[s][0] += tm; a[s][1] += tv; a[s][2] += n

    if engine == "modkit":
        for k in keys:
            ids = set(reads.read_id[(reads.label == k[0]) & (reads.quintile == k[1])])
            arr = P.pileup_modkit(bam_paths, thresholds, filt, label_by, modkit=modkit,
                                  threads=threads, read_ids=ids)
            for chrom, mv2 in _strand_sum(arr).items():
                consume(k, chrom, mv2)
    else:
        qmap = dict(zip(zip(reads.label, reads.read_id), reads.quintile))
        lens = P.chrom_lengths(bam_paths[0])
        for chrom in sorted(reads.chrom.unique()):
            cur = {k: np.zeros((2, lens[chrom]), np.int32) for k in keys}
            for read, lab, _ in iter_reads(bam_paths, (chrom, 0, lens[chrom]), filt, label_by):
                q = qmap.get((lab, read.query_name), 0)
                if not q:
                    continue
                f = read_to_fiber(read, thresholds.get(lab), lab)
                called = f.call >= 0
                pos = f.pos[called]
                a = cur[(lab, q)]
                np.add.at(a[1], pos, 1)
                np.add.at(a[0], pos[f.call[called] == 1], 1)
            for k in keys:
                consume(k, chrom, cur[k])
            if progress:
                print(f"[quintiles] pass 2: {chrom} done", file=sys.stderr)
    for (lab, q) in keys:
        sel = reads[(reads.label == lab) & (reads.quintile == q)]
        a = acc[(lab, q)]
        summ.append(dict(label=lab, quintile=q, n_reads=len(sel), rate_lo=sel.rate.min(),
                         rate_hi=sel.rate.max(), rate_median=sel.rate.median(),
                         genome_n_mod=a["gm"], genome_n_valid=a["gv"],
                         genome_m6a_rate=a["gm"] / a["gv"] if a["gv"] else np.nan))
        for s in site_specs:
            h = halves[s]
            tm, tv, n = a[s]
            prof_rows.append(pd.DataFrame({
                "label": lab, "quintile": q, "site_class": refs.site_label(s),
                "pos": np.arange(-h, h + 1), "n_mod": tm, "n_valid": tv, "n_sites": n}))
    prof = pd.concat(prof_rows, ignore_index=True)
    with np.errstate(all="ignore"):
        prof["m6a_rate"] = prof.n_mod / prof.n_valid
    return reads, prof, pd.DataFrame(summ)


# --------------------------------------------------------------------------
# footprint scores (from footprint_metaprofiles.py)
# --------------------------------------------------------------------------
def _mean_abs(y, x, lo, hi):
    return np.nanmean(y[(np.abs(x) >= lo) & (np.abs(x) <= hi)])


def _mean_between(y, x, lo, hi):
    return np.nanmean(y[(x >= lo) & (x <= hi)])


def nrl_power(y, x):
    d = y - pd.Series(y).rolling(401, center=True, min_periods=1).mean().values
    d = np.nan_to_num(d)
    power = np.abs(np.fft.rfft(d)) ** 2
    freqs = np.fft.rfftfreq(len(d), d=1.0)
    with np.errstate(divide="ignore"):
        periods = 1 / freqs
    band = (periods >= 150) & (periods <= 185)
    keep = (periods >= 20) & (periods <= 600)
    return power[band].sum() / power[keep].sum()


def scores(prof: pd.DataFrame, summ: pd.DataFrame, smooth_bp: int = 11) -> pd.DataFrame:
    out = []
    for (lab, q), g in prof.groupby(["label", "quintile"]):
        gr = summ[(summ.label == lab) & (summ.quintile == q)].genome_m6a_rate.iloc[0]
        s = dict(label=lab, quintile=q, genome_m6a_rate=gr)
        for cls, d in g.groupby("site_class"):
            y = pd.Series(d.m6a_rate.values / gr).rolling(smooth_bp, center=True, min_periods=1).mean().values
            x = d.pos.values
            if cls in ("dyad", "nuc"):
                s["dyad_linker_over_dyad"] = _mean_abs(y, x, 80, 100) / _mean_abs(y, x, 0, 20)
                s["dyad_nrl_power"] = nrl_power(y, x)
            elif cls == "tss":
                s["tss_ndr_over_body"] = _mean_between(y, x, -200, 0) / _mean_between(y, x, 500, 1500)
            else:
                s[f"{cls}_flank_over_center"] = _mean_abs(y, x, 30, 100) / _mean_abs(y, x, 0, 10)
        out.append(s)
    return pd.DataFrame(out)
