"""Per-sample m6A call thresholds ("10th-percentile confidence", modkit style).

Algorithm (identical to `modkit sample-probs` / the roboNhat per-barcode TSV):

* every m6A target base in the stored read (A on + reads, T on - reads,
  soft-clips included) contributes one call;
* p_mod = (ML + 0.5) / 256, confidence = max(p_mod, 1 - p_mod);
* bases absent from MM count as canonical with confidence 1.0 when the MM mode
  is '.', and are skipped when it is '?';
* threshold = 10th percentile of the pooled confidence distribution, which
  always lands on a bin value (k + 0.5) / 256;
* ml_meth = round(thr * 256 - 0.5), ml_canon = 255 - ml_meth.

Two engines:
  pysam  - pure Python (default, no extra binaries)
  modkit - runs `modkit sample-probs` per label (needs one label per BAM file)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from .bam import (CallThreshold, ReadFilter, Thresholds, iter_reads, m6a_ml,
                  open_bam, read_label, target_mask)

N_BINS = 257  # ML 0..255 + one bin for unlisted ('.' mode) A's
UNLISTED = 256


def modkit_conf() -> np.ndarray:
    p = (np.arange(256) + 0.5) / 256.0
    return np.append(np.maximum(p, 1 - p), 1.0)


def percentile_thr(h: np.ndarray, q: float = 0.10) -> float:
    """Lowest confidence c with P(conf <= c) >= q."""
    c = modkit_conf()
    o = np.argsort(c, kind="stable")
    cum = np.cumsum(h[o]) / h.sum()
    return float(c[o][np.searchsorted(cum, q)])


def add_read(h: np.ndarray, read) -> int:
    """Add one read's target-base ML values to histogram h. Returns #bases."""
    ml, mode = m6a_ml(read)
    v = ml[target_mask(read)]
    listed = v[v >= 0]
    h[:256] += np.bincount(listed, minlength=256)[:256]
    if mode != "?":
        h[UNLISTED] += int((v < 0).sum())
    return int(v.size)


def histograms(bam_paths: Sequence[str], filt: Optional[ReadFilter] = None,
               label_by: str = "auto", max_reads_per_label: Optional[int] = 50000,
               seed: int = 0) -> Dict[str, dict]:
    """Per-label 257-bin ML histograms.

    max_reads_per_label: stop adding reads to a label after this many (reads
    are taken in file order; set None to use every read, as the lab TSV did).
    """
    filt = filt if filt is not None else ReadFilter()
    out: Dict[str, dict] = defaultdict(lambda: {"h": np.zeros(N_BINS, np.int64), "n_reads": 0})
    for read, label, _ in iter_reads(bam_paths, None, filt, label_by):
        d = out[label]
        if max_reads_per_label is not None and d["n_reads"] >= max_reads_per_label:
            continue
        add_read(d["h"], read)
        d["n_reads"] += 1
    return dict(out)


def estimate(bam_paths: Sequence[str], filt: Optional[ReadFilter] = None,
             label_by: str = "auto", percentile: float = 0.10,
             max_reads_per_label: Optional[int] = 50000, engine: str = "pysam",
             modkit: Optional[str] = None, threads: int = 4) -> pd.DataFrame:
    """Estimate per-label thresholds. Returns a DataFrame (see write_tsv)."""
    if engine == "modkit":
        return estimate_modkit(bam_paths, filt, label_by, percentile,
                               max_reads_per_label, modkit, threads)
    rows = []
    for label, d in sorted(histograms(bam_paths, filt, label_by, max_reads_per_label).items()):
        h = d["h"]
        if h.sum() == 0:
            continue
        thr = percentile_thr(h, percentile)
        ct = CallThreshold.from_conf(thr)
        rows.append(dict(label=label, thr=thr, ml_meth=ct.ml_meth, ml_canon=ct.ml_canon,
                         n_reads=d["n_reads"], n_calls=int(h.sum()),
                         frac_unlisted=float(h[UNLISTED] / h.sum()),
                         percentile=percentile, engine="pysam"))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# modkit engine
# --------------------------------------------------------------------------
def find_modkit(modkit: Optional[str] = None) -> str:
    cand = modkit or os.environ.get("FIBERNHAT_MODKIT") or shutil.which("modkit")
    if not cand or not os.path.exists(cand):
        raise FileNotFoundError(
            "modkit not found: install `ont-modkit` (bioconda), put it on PATH, "
            "set FIBERNHAT_MODKIT, or set tools.modkit in config/site.yaml")
    return cand


def write_filtered_bam(bam_path: str, out_path: str, filt: ReadFilter,
                       label_by: str = "auto", only_label: Optional[str] = None,
                       read_ids: Optional[set] = None) -> int:
    """Copy reads passing `filt` (and optionally one label / a read-id set) to
    a new indexed BAM. Used to feed modkit exactly the same reads as pysam."""
    n = 0
    with open_bam(bam_path) as bam, \
            __import__("pysam").AlignmentFile(out_path, "wb", template=bam) as out:
        for read in bam.fetch(until_eof=True):
            if filt.reason(read) is not None:
                continue
            if only_label is not None and read_label(read, bam_path, label_by) != only_label:
                continue
            if read_ids is not None and read.query_name not in read_ids:
                continue
            out.write(read)
            n += 1
    import pysam
    pysam.sort("-o", out_path + ".sorted.bam", out_path)
    os.replace(out_path + ".sorted.bam", out_path)
    pysam.index(out_path)
    return n


def _labels_in(bam_path: str, filt: ReadFilter, label_by: str) -> List[str]:
    labs = set()
    for _, lab, _ in iter_reads([bam_path], None, filt, label_by):
        labs.add(lab)
    return sorted(labs)


def parse_modkit_thresholds(path: str, percentile: float = 0.10) -> float:
    df = pd.read_csv(path, sep=r"\s+")
    df.columns = [c.lower() for c in df.columns]
    df = df[df["base"] == "A"]
    row = df.iloc[(df["percentile"] - percentile).abs().argsort().iloc[0]]
    return float(row["threshold"])


def estimate_modkit(bam_paths, filt, label_by, percentile, max_reads_per_label,
                    modkit, threads) -> pd.DataFrame:
    mk = find_modkit(modkit)
    filt = filt if filt is not None else ReadFilter()
    rows = []
    with tempfile.TemporaryDirectory(prefix="fibernhat_mk_") as tmp:
        for p in bam_paths:
            for label in _labels_in(p, filt, label_by):
                fb = os.path.join(tmp, f"{label}.bam")
                n = write_filtered_bam(p, fb, filt, label_by, only_label=label)
                od = os.path.join(tmp, f"sp_{label}")
                cmd = [mk, "sample-probs", fb, "-t", str(threads), "-o", od,
                       "--percentiles", str(percentile), "--force", "--suppress-progress"]
                if max_reads_per_label is None or n <= max_reads_per_label:
                    cmd += ["--no-sampling"]
                else:
                    cmd += ["-n", str(max_reads_per_label), "--seed", "0"]
                subprocess.run(cmd, check=True, capture_output=True)
                thr = parse_modkit_thresholds(os.path.join(od, "thresholds.tsv"), percentile)
                ct = CallThreshold.from_conf(thr)
                rows.append(dict(label=label, thr=thr, ml_meth=ct.ml_meth, ml_canon=ct.ml_canon,
                                 n_reads=n, n_calls=np.nan, frac_unlisted=np.nan,
                                 percentile=percentile, engine="modkit"))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# TSV io
# --------------------------------------------------------------------------
def write_tsv(df: pd.DataFrame, path: str, note: str = "") -> None:
    with open(path, "w") as fh:
        fh.write("# fibernhat per-label m6A call thresholds (10th-percentile confidence)\n")
        fh.write("# call: meth if ML >= ml_meth; canonical if ML <= ml_canon or A unlisted ('.'); else no-call\n")
        if note:
            fh.write(f"# {note}\n")
        df.to_csv(fh, sep="\t", index=False, float_format="%.7f")


def read_tsv(path: str) -> Thresholds:
    """Load a thresholds TSV (fibernhat format, or the lab roboNhat format
    with 'barcode'/'rg' columns -- both the barcode and RG become keys)."""
    df = pd.read_csv(path, sep="\t", comment="#")
    by = {}
    for _, r in df.iterrows():
        ct = CallThreshold(int(r["ml_meth"]), int(r["ml_canon"]))
        for key in ("label", "barcode", "rg"):
            if key in df.columns and isinstance(r[key], str):
                by[r[key]] = ct
    return Thresholds(by, None, source=os.path.basename(path))


def resolve(thresholds_tsv: Optional[str], threshold_ml: Optional[int],
            bam_paths: Sequence[str], filt: Optional[ReadFilter], label_by: str,
            max_reads_per_label: Optional[int] = 20000, verbose: bool = True) -> Thresholds:
    """Shared CLI logic: fixed ML > TSV > estimate on the fly from the BAMs."""
    if threshold_ml is not None:
        return Thresholds.fixed(threshold_ml)
    if thresholds_tsv:
        return read_tsv(thresholds_tsv)
    df = estimate(bam_paths, filt, label_by, max_reads_per_label=max_reads_per_label)
    if verbose:
        import sys
        print("[fibernhat] estimated per-label thresholds:\n" +
              df[["label", "thr", "ml_meth", "ml_canon", "n_reads"]].to_string(index=False),
              file=sys.stderr)
    return Thresholds({r.label: CallThreshold(int(r.ml_meth), int(r.ml_canon))
                       for r in df.itertuples()}, None, source="estimated")
