"""Per-read QC table straight from BAM (aligned or unaligned dorado output).

One row per primary record (unmapped reads included; secondary and
supplementary records skipped so each read counts once):
    read_id, label, mapped, chrom, mapq, read_len, ref_span, mean_q, has_mm,
    fail_reason (first failing step of the read filter, '' = passes)
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..io.bam import (CASCADE_STEPS, ReadFilter, open_bam, read_label, read_length,
                      read_mean_q)


def read_table(bam_paths: Sequence[str], label_by: str = "auto",
               filt: Optional[ReadFilter] = None, max_reads: Optional[int] = None) -> pd.DataFrame:
    filt = filt if filt is not None else ReadFilter()
    rows = []
    for p in bam_paths:
        n = 0
        with open_bam(p) as bam:
            for r in bam.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary:
                    continue
                q = read_mean_q(r)
                rows.append((r.query_name, read_label(r, p, label_by), not r.is_unmapped,
                             r.reference_name if not r.is_unmapped else "*",
                             r.mapping_quality, read_length(r),
                             (r.reference_end - r.reference_start) if not r.is_unmapped else 0,
                             np.nan if q is None else q,
                             r.has_tag("MM") or r.has_tag("Mm"), filt.reason(r) or ""))
                n += 1
                if max_reads and n >= max_reads:
                    break
    return pd.DataFrame(rows, columns=["read_id", "label", "mapped", "chrom", "mapq", "read_len",
                                       "ref_span", "mean_q", "has_mm", "fail_reason"])


def n50(lengths) -> int:
    a = np.sort(np.asarray(lengths))[::-1]
    if a.size == 0:
        return 0
    c = np.cumsum(a)
    return int(a[np.searchsorted(c, c[-1] / 2)])


def summary(df: pd.DataFrame, genome_size: int) -> pd.DataFrame:
    out = []
    for lab, g in df.groupby("label", sort=False):
        L = g.read_len.values
        passed = g[g.fail_reason == ""]
        out.append(dict(
            label=lab, n_reads=len(g), total_bases=int(L.sum()), n50=n50(L),
            median_len=float(np.median(L)) if len(L) else 0, mean_len=float(L.mean()) if len(L) else 0,
            max_len=int(L.max()) if len(L) else 0,
            median_q=float(np.nanmedian(g.mean_q)), mean_q=float(np.nanmean(g.mean_q)),
            frac_q10=float((g.mean_q >= 10).mean()), frac_mapped=float(g.mapped.mean()),
            frac_with_mm=float(g.has_mm.mean()), coverage_x=float(L.sum() / genome_size),
            n_pass_filter=len(passed), pass_bases=int(passed.read_len.sum()),
            pass_coverage_x=float(passed.ref_span.sum() / genome_size)))
    return pd.DataFrame(out)


def cascade(df: pd.DataFrame) -> pd.DataFrame:
    """Reads (and bases) removed at each filter step, per label."""
    rows = []
    for lab, g in df.groupby("label", sort=False):
        remaining, rb = len(g), int(g.read_len.sum())
        rows.append(dict(label=lab, step="all_primary", removed=0, removed_bases=0,
                         remaining=remaining, remaining_bases=rb))
        for step in CASCADE_STEPS:
            s = g[g.fail_reason == step]
            if step == "not_in_keep_list" and s.empty:
                continue
            remaining -= len(s); rb -= int(s.read_len.sum())
            rows.append(dict(label=lab, step=step, removed=len(s), removed_bases=int(s.read_len.sum()),
                             remaining=remaining, remaining_bases=rb))
    return pd.DataFrame(rows)
