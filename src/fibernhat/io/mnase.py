"""Paired-end MNase-seq fragments from BAM (vendored from DMS_seq read_bam.py).

One fragment per read pair: keep the read whose mate is reverse (== samtools
-f 32), length = template_length, midpoint = pos + 1 + length/2 (1-based).
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd
import pysam

from .refs import ROMAN


def _fetch(bam: pysam.AlignmentFile, chrom: str, start: int, end: int):
    names = set(bam.references)
    for c in (chrom, chrom.replace("chr", ""), ):
        if c in names:
            return bam.fetch(c, max(0, start), end)
    # numeric chrom names (1..16)
    roman = chrom.replace("chr", "")
    if roman in ROMAN and str(ROMAN.index(roman) + 1) in names:
        return bam.fetch(str(ROMAN.index(roman) + 1), max(0, start), end)
    return iter(())


def fragments(bam_path: str, chrom: str, start: int, end: int,
              max_len: int = 300, pad: int = 300) -> pd.DataFrame:
    """Fragments whose midpoint falls in [start, end) (0-based); columns
    mid (1-based genomic bp), length."""
    mids, lens = [], []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in _fetch(bam, chrom, start - pad, end + pad):
            if not r.mate_is_reverse or r.is_reverse:
                continue
            L = r.template_length
            if L <= 0 or L > max_len:
                continue
            mid = r.reference_start + 1 + L / 2
            if start <= mid - 1 < end:
                mids.append(mid)
                lens.append(L)
    return pd.DataFrame({"mid": np.asarray(mids, float), "length": np.asarray(lens, int)})


def fragments_around_sites(bam_path: str, sites: pd.DataFrame, half: int,
                           max_len: int = 250, max_points: Optional[int] = 10000,
                           seed: int = 0) -> pd.DataFrame:
    """Pool fragments around many sites (V-plot); x = oriented distance."""
    rel, ln = [], []
    for c, g in sites.groupby("chrom"):
        for center, sign in zip(g.center.values, g.sign.values):
            f = fragments(bam_path, c, center - half, center + half + 1, max_len, pad=max_len)
            if len(f):
                rel.append((f.mid.values - 1 - center) * sign)
                ln.append(f.length.values)
    if not rel:
        return pd.DataFrame({"rel": [], "length": []})
    df = pd.DataFrame({"rel": np.concatenate(rel), "length": np.concatenate(ln)})
    if max_points and len(df) > max_points:
        df = df.sample(max_points, random_state=seed)
    return df
