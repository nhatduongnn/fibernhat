"""Per-position m6A pileup (bedMethyl-like counts) from BAM.

pysam engine: for every read passing the filter, every aligned target base
(A on + reads, T on - reads) with a call contributes Nvalid += 1 at its
reference position on that strand; methylated calls also add Nmod += 1.
No-calls (inside the per-label confidence band) are skipped, exactly like
`modkit pileup --filter-threshold A:<thr>`; with a fixed ML>=128 binary call
it equals `modkit pileup --no-filtering`.

modkit engine: writes a filtered BAM per label and runs `modkit pileup`.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .io.bam import (ReadFilter, Thresholds, iter_reads, open_bam, read_to_fiber)

# arrays[chrom] = int32[4, L]: rows = (mod +, valid +, mod -, valid -)
Arrays = Dict[str, np.ndarray]


def empty_arrays(chrom_lens: Dict[str, int]) -> Arrays:
    return {c: np.zeros((4, L), dtype=np.int32) for c, L in chrom_lens.items()}


def add_fiber(arrays: Arrays, f) -> None:
    a = arrays.get(f.chrom)
    if a is None:
        return
    called = f.call >= 0
    pos = f.pos[called]
    meth = f.call[called] == 1
    o = 0 if f.strand == "+" else 2
    np.add.at(a[o + 1], pos, 1)
    np.add.at(a[o], pos[meth], 1)


def chrom_lengths(bam_path: str) -> Dict[str, int]:
    with open_bam(bam_path) as b:
        return dict(zip(b.references, b.lengths))


def pileup_pysam(bam_paths: Sequence[str], thresholds: Thresholds,
                 filt: Optional[ReadFilter] = None, label_by: str = "auto",
                 region: Optional[Tuple[str, int, int]] = None,
                 read_ids: Optional[set] = None) -> Arrays:
    filt = filt if filt is not None else ReadFilter()
    lens = chrom_lengths(bam_paths[0])
    if region is not None:
        lens = {region[0]: lens[region[0]]}
    arrays = empty_arrays(lens)
    for read, label, _ in iter_reads(bam_paths, region, filt, label_by):
        if read_ids is not None and read.query_name not in read_ids:
            continue
        add_fiber(arrays, read_to_fiber(read, thresholds.get(label), label))
    return arrays


def arrays_to_frame(arrays: Arrays) -> pd.DataFrame:
    rows = []
    for c, a in arrays.items():
        for strand, o in (("+", 0), ("-", 2)):
            idx = np.nonzero(a[o + 1])[0]
            if idx.size:
                rows.append(pd.DataFrame({"chrom": c, "pos": idx, "strand": strand,
                                          "n_mod": a[o, idx], "n_valid": a[o + 1, idx]}))
    if not rows:
        return pd.DataFrame(columns=["chrom", "pos", "strand", "n_mod", "n_valid"])
    return pd.concat(rows, ignore_index=True)


def write_bedmethyl_like(arrays: Arrays, path: str) -> None:
    df = arrays_to_frame(arrays)
    df["end"] = df.pos + 1
    df["frac"] = 100 * df.n_mod / df.n_valid
    df[["chrom", "pos", "end", "strand", "n_valid", "n_mod", "frac"]].to_csv(
        path, sep="\t", index=False, header=False, float_format="%.2f")


# --------------------------------------------------------------------------
# modkit engine
# --------------------------------------------------------------------------
def read_modkit_bedmethyl(path: str, chrom_lens: Dict[str, int]) -> Arrays:
    cols = [0, 1, 3, 5, 9, 11]
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=cols,
                     names=["chrom", "start", "code", "strand", "n_valid", "n_mod"])
    df = df[df.code.astype(str).isin(["a", "21891"])]
    arrays = empty_arrays({c: chrom_lens[c] for c in df.chrom.unique() if c in chrom_lens})
    for (c, s), g in df.groupby(["chrom", "strand"]):
        if c not in arrays:
            continue
        o = 0 if s == "+" else 2
        np.add.at(arrays[c][o], g.start.values, g.n_mod.values)
        np.add.at(arrays[c][o + 1], g.start.values, g.n_valid.values)
    return arrays


def modkit_pileup_bam(modkit: str, bam: str, out_bed: str, conf_thr: Optional[float],
                      threads: int = 4, region: Optional[str] = None) -> None:
    """conf_thr=None -> --no-filtering (argmax call, equals fixed ML>=128)."""
    cmd = [modkit, "pileup", bam, out_bed, "-t", str(threads), "--suppress-progress"]
    if conf_thr is None:
        cmd += ["--no-filtering"]
    else:
        cmd += ["--filter-threshold", f"A:{conf_thr}"]
    if region:
        cmd += ["--region", region]
    subprocess.run(cmd, check=True, capture_output=True)


def pileup_modkit(bam_paths: Sequence[str], thresholds: Thresholds,
                  filt: Optional[ReadFilter] = None, label_by: str = "auto",
                  modkit: Optional[str] = None, threads: int = 4,
                  read_ids: Optional[set] = None) -> Arrays:
    """Same output as pileup_pysam, computed by modkit. Each label is written
    to its own filtered BAM so per-label thresholds can be applied.

    Fixed-threshold mode is only exact for ML>=128 (modkit's argmax call)."""
    from .io.thresholds import find_modkit, write_filtered_bam, _labels_in
    mk = find_modkit(modkit)
    filt = filt if filt is not None else ReadFilter()
    lens = chrom_lengths(bam_paths[0])
    total = empty_arrays(lens)
    with tempfile.TemporaryDirectory(prefix="fibernhat_pu_") as tmp:
        for i, p in enumerate(bam_paths):
            for label in _labels_in(p, filt, label_by):
                ct = thresholds.get(label)
                fixed = ct.ml_canon == ct.ml_meth - 1
                if fixed and ct.ml_meth != 128:
                    raise ValueError("modkit engine supports fixed threshold ML>=128 only "
                                     "(modkit argmax); use per-label thresholds or pysam")
                fb = os.path.join(tmp, f"{i}_{label}.bam")
                write_filtered_bam(p, fb, filt, label_by, only_label=label, read_ids=read_ids)
                bed = fb + ".bed"
                modkit_pileup_bam(mk, fb, bed, None if fixed else ct.conf, threads)
                for c, a in read_modkit_bedmethyl(bed, lens).items():
                    total[c] += a
    return total


def pileup(bam_paths, thresholds, filt=None, label_by="auto", engine="pysam", **kw) -> Arrays:
    if engine == "modkit":
        return pileup_modkit(bam_paths, thresholds, filt, label_by,
                             modkit=kw.get("modkit"), threads=kw.get("threads", 4),
                             read_ids=kw.get("read_ids"))
    return pileup_pysam(bam_paths, thresholds, filt, label_by,
                        region=kw.get("region"), read_ids=kw.get("read_ids"))


def compare(a: Arrays, b: Arrays) -> pd.DataFrame:
    """Positions where two pileups disagree (empty frame = identical)."""
    fa, fb = arrays_to_frame(a), arrays_to_frame(b)
    m = fa.merge(fb, on=["chrom", "pos", "strand"], how="outer", suffixes=("_a", "_b")).fillna(0)
    return m[(m.n_mod_a != m.n_mod_b) | (m.n_valid_a != m.n_valid_b)]
