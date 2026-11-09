"""Reference files shipped in ref/ and site-set loaders.

All site loaders return a DataFrame with columns
    chrom (UCSC style, e.g. chrIV), center (0-based int), sign (+1 / -1), name
`sign = -1` means the site is on the Crick strand and windows are flipped so
the site's 5'->3' direction always points right.
"""
from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
REF_DIR = Path(os.environ.get("FIBERNHAT_REF", REPO_ROOT / "ref"))

FASTA = REF_DIR / "sacCer3" / "custom_sacCer3_lambda.fasta"
CHROM_SIZES = REF_DIR / "sacCer3" / "sacCer3_lambda.chrom.sizes"
SGD_GFF = REF_DIR / "sacCer3" / "saccharomyces_cerevisiae_R64-1-1_20110208_no_fasta.gff"
TSS_GFF = REF_DIR / "sacCer3" / "Park_2014_TSS_V64.gff"
PAS_GFF = REF_DIR / "sacCer3" / "Park_2014_PAS_V64.gff"
SITES = REF_DIR / "sites"

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
         "XI", "XII", "XIII", "XIV", "XV", "XVI"]
ARABIC_TO_ROMAN = {f"chr{i + 1}": f"chr{r}" for i, r in enumerate(ROMAN)}
NUCLEAR = [f"chr{r}" for r in ROMAN]


def norm_chrom(c: str) -> str:
    """'chr9' / '9' / 'IX' / 'chrIX' -> 'chrIX'; 'chrmt'/'chrMito' -> 'chrM'."""
    c = str(c)
    if c in ARABIC_TO_ROMAN:
        return ARABIC_TO_ROMAN[c]
    if c.isdigit() and f"chr{c}" in ARABIC_TO_ROMAN:
        return ARABIC_TO_ROMAN[f"chr{c}"]
    if c in ROMAN:
        return f"chr{c}"
    if c.lower() in ("chrmt", "chrmito", "mt", "mito", "chrm"):
        return "chrM"
    return c


@lru_cache(maxsize=4)
def chrom_sizes(path: str = str(CHROM_SIZES)) -> Dict[str, int]:
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["chrom", "len"])
    return dict(zip(df.chrom, df.len.astype(int)))


def nuclear_genome_size(path: str = str(CHROM_SIZES)) -> int:
    return int(sum(v for k, v in chrom_sizes(path).items() if k in NUCLEAR))


def _frame(chrom, center, sign, name=None) -> pd.DataFrame:
    df = pd.DataFrame({"chrom": [norm_chrom(c) for c in chrom],
                       "center": np.asarray(center).astype(int),
                       "sign": np.asarray(sign).astype(int)})
    df["name"] = list(name) if name is not None else [f"{c}:{p + 1}" for c, p in zip(df.chrom, df.center)]
    return df.drop_duplicates(["chrom", "center", "sign"]).reset_index(drop=True)


def load_dyads(path: str = str(SITES / "Brogaard_nuc_positions.sacCer3.top2000.tsv")) -> pd.DataFrame:
    """Brogaard 2012 top-2000 nucleosome dyads (1-based positions)."""
    d = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1], names=["chrom", "pos"])
    return _frame(d.chrom, d.pos - 1, np.ones(len(d)))


def load_tss(path: str = str(TSS_GFF)) -> pd.DataFrame:
    """Park 2014 TSSs, oriented by gene strand (GFF col 4 is 1-based)."""
    g = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 6],
                    names=["chrom", "name", "pos", "strand"])
    g = g[g.pos.astype(str).str.isdigit()]
    return _frame(g.chrom, g.pos.astype(int) - 1, np.where(g.strand == "-", -1, 1),
                  g.name.str.replace("_TSS", "", regex=False))


def load_fimo(tf: str, sites_dir: str = str(SITES)) -> pd.DataFrame:
    """Rossi 2021 ChIP-exo summits with a FIMO motif (abf1/reb1 shipped),
    centered on the motif midpoint, oriented by motif strand."""
    b = pd.read_csv(os.path.join(sites_dir, f"{tf.lower()}_fimo.bed"), sep="\t")
    b = b[b.has_motif.astype(str) == "True"]
    b = b.astype({"motif_start": int, "motif_stop": int})
    b = b.drop_duplicates(["motif_chrom", "motif_start", "motif_stop"])
    return _frame(b.motif_chrom, (b.motif_start + b.motif_stop) // 2,
                  np.where(b.motif_strand == "-", -1, 1), b.name)


def load_chereji(path: str = str(SITES / "Chereji_2018_+1_-1_nucs.bed")) -> pd.DataFrame:
    b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chrom", "start"])
    return _frame(b.chrom, b.start, np.ones(len(b)))


def load_rossi_all(path: str = str(SITES / "rossi_peak_w_strand_conformed_to_PWM_all_TFs_peakVal.bed"),
                   tf: str = None) -> pd.DataFrame:
    b = pd.read_csv(path, sep="\t")
    if tf:
        b = b[b.TF.str.lower().str.startswith(tf.lower())]
    return _frame(b.chr, (b.start + b.end) // 2, np.where(b.strand == "-", -1, 1), b.TF)


def load_generic(path: str) -> pd.DataFrame:
    """TSV/BED with a header containing chrom + (pos|dyad|center|start[,end])
    + optional strand. Headerless BED (chrom start end [name score strand]) also
    works. Positions in a 'pos'/'dyad' column are 1-based; BED starts 0-based."""
    with open(path) as fh:
        first = fh.readline().lstrip("#").strip().split("\t")
    has_header = not first[1].strip().isdigit() if len(first) > 1 else False
    if has_header:
        df = pd.read_csv(path, sep="\t", comment=None)
        df.columns = [c.lstrip("#").lower() for c in df.columns]
        chrom = df[[c for c in df.columns if c in ("chrom", "chr", "chromosome")][0]]
        if "center" in df.columns:
            center = df.center
        elif "pos" in df.columns or "dyad" in df.columns:
            center = (df["pos"] if "pos" in df.columns else df["dyad"]) - 1
        else:
            st = df["start"]
            center = (st + df["end"]) // 2 if "end" in df.columns else st
        strand = df["strand"] if "strand" in df.columns else pd.Series(["+"] * len(df))
    else:
        df = pd.read_csv(path, sep="\t", header=None)
        chrom, st = df[0], df[1]
        center = (st + df[2] - 1) // 2 if df.shape[1] > 2 else st
        strand = df[5] if df.shape[1] > 5 else pd.Series(["+"] * len(df))
    return _frame(chrom, center, np.where(strand.astype(str) == "-", -1, 1))


SITE_SETS = {
    "dyad": load_dyads,
    "nuc": load_dyads,
    "tss": load_tss,
    "abf1": lambda: load_fimo("abf1"),
    "reb1": lambda: load_fimo("reb1"),
    "chereji": load_chereji,
    "tf": load_rossi_all,
}

# default half-windows (bp) per site class, from footprint_metaprofiles.py
DEFAULT_HALF = {"dyad": 1000, "nuc": 1000, "tss": 1500, "abf1": 300, "reb1": 300,
                "chereji": 1000, "tf": 300}


def load_sites(spec: str) -> pd.DataFrame:
    """'dyad' | 'tss' | 'abf1' | 'reb1' | 'chereji' | 'tf' | 'tf:Reb1' | path."""
    if spec in SITE_SETS:
        return SITE_SETS[spec]()
    if spec.startswith("tf:"):
        return load_rossi_all(tf=spec.split(":", 1)[1])
    if os.path.exists(spec):
        return load_generic(spec)
    raise ValueError(f"unknown site set {spec!r}; use one of {sorted(SITE_SETS)} or a file path")


def site_label(spec: str) -> str:
    return Path(spec).stem if os.path.exists(spec) else spec
