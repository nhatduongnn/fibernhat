"""Read nanopore fiber-seq (m6A) BAMs directly with pysam.

This is the single place where MM/ML tags are parsed, reads are filtered and
m6A calls are made. Every skill (QC, fiberplot, aggregate, quintiles) goes
through these functions so a threshold or filter change applies everywhere.

Conventions
-----------
* Reference coordinates are 0-based, half-open (pysam convention). Region
  strings ``chrI:60500-63500`` are parsed as 1-based inclusive (samtools
  convention) and converted.
* m6A is the ``A+a`` (or ChEBI 21891) MM group. For reverse-strand reads the
  stored SEQ shows the methylated template A as ``T``; calls are reported at the
  reference position of that base (Crick strand = ref ``T``).
* ML byte -> probability uses the modkit convention ``p = (ML + 0.5) / 256``.
* Call states: 1 = methylated, 0 = canonical, -1 = no-call.
  ``ML >= ml_meth`` -> 1; ``ML <= ml_canon`` -> 0; A absent from MM -> 0 when
  the MM mode is ``.`` (implicit canonical), -1 when the mode is ``?``.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pysam

M6A_CODES = ("a", 21891)

# rDNA locus on chrXII (sacCer3), excluded by default (huge copy number).
RDNA = ("chrXII", 451417, 468929)

_BARCODE_RE = re.compile(r"(barcode\d+|unclassified)", re.IGNORECASE)


# --------------------------------------------------------------------------
# regions
# --------------------------------------------------------------------------
def parse_region(region: str) -> Tuple[str, int, int]:
    """'chrI:60,500-63,500' -> ('chrI', 60499, 63500)  (0-based half-open)."""
    region = region.replace(",", "")
    if ":" not in region:
        raise ValueError(f"region must look like chrI:1000-2000, got {region!r}")
    chrom, span = region.rsplit(":", 1)
    s, e = span.split("-")
    s, e = int(s), int(e)
    if e <= s - 1:
        raise ValueError(f"region end before start: {region}")
    return chrom, s - 1, e


# --------------------------------------------------------------------------
# read filters
# --------------------------------------------------------------------------
@dataclass
class ReadFilter:
    """Default = the lab 'filtnordna' cascade (roboNhat, 2026-11)."""
    min_mapq: int = 20
    min_span: int = 500             # reference span (reference_end - reference_start)
    min_qs: Optional[float] = 10.0  # dorado mean-Q tag 'qs'; None disables
    exclude_chroms: Tuple[str, ...] = ("chrM", "chrLambda")
    exclude_regions: Tuple[Tuple[str, int, int], ...] = (RDNA,)
    require_mm: bool = True
    keep_read_ids: Optional[set] = None

    @classmethod
    def none(cls) -> "ReadFilter":
        """Primary mapped reads only, nothing else (modkit-like)."""
        return cls(min_mapq=0, min_span=0, min_qs=None, exclude_chroms=(),
                   exclude_regions=(), require_mm=False)

    @classmethod
    def from_dict(cls, d: Optional[dict]) -> "ReadFilter":
        if not d:
            return cls()
        d = dict(d)
        if "exclude_regions" in d and d["exclude_regions"] is not None:
            regs = []
            for r in d["exclude_regions"]:
                regs.append(parse_region(r) if isinstance(r, str) else tuple(r))
            d["exclude_regions"] = tuple(regs)
        if "exclude_chroms" in d and d["exclude_chroms"] is not None:
            d["exclude_chroms"] = tuple(d["exclude_chroms"])
        return cls(**d)

    def reason(self, read: pysam.AlignedSegment) -> Optional[str]:
        """Return why a read fails (first failing step), or None if it passes.

        Order matches the cascade reported by `fibernhat qc`.
        """
        if read.is_unmapped:
            return "unmapped"
        if read.is_secondary or read.is_supplementary:
            return "secondary_supplementary"
        if read.is_qcfail or read.is_duplicate:
            return "qcfail_duplicate"
        if self.keep_read_ids is not None and read.query_name not in self.keep_read_ids:
            return "not_in_keep_list"
        if self.require_mm and not (read.has_tag("MM") or read.has_tag("Mm")):
            return "no_mm_tag"
        chrom = read.reference_name
        if chrom in self.exclude_chroms:
            return "excluded_chrom"
        for c, s, e in self.exclude_regions:
            if chrom == c and read.reference_start < e and read.reference_end > s:
                return "excluded_region"
        if read.mapping_quality < self.min_mapq:
            return "low_mapq"
        if (read.reference_end - read.reference_start) < self.min_span:
            return "short_span"
        if self.min_qs is not None:
            qs = read_mean_q(read)
            if qs is not None and qs < self.min_qs:
                return "low_qs"
        return None

    def passes(self, read: pysam.AlignedSegment) -> bool:
        return self.reason(read) is None


CASCADE_STEPS = ("unmapped", "secondary_supplementary", "qcfail_duplicate",
                 "not_in_keep_list", "no_mm_tag", "excluded_chrom",
                 "excluded_region", "low_mapq", "short_span", "low_qs")


# --------------------------------------------------------------------------
# per-read helpers
# --------------------------------------------------------------------------
def read_mean_q(read: pysam.AlignedSegment) -> Optional[float]:
    """Dorado 'qs' tag if present, else mean Q from base qualities.

    The fallback averages error probabilities (like dorado) over all bases;
    dorado itself skips the first 60 bases, so values can differ slightly.
    """
    if read.has_tag("qs"):
        return float(read.get_tag("qs"))
    q = read.query_qualities
    if q is None or len(q) == 0:
        return None
    q = np.asarray(q, dtype=float)
    return float(-10.0 * np.log10(np.mean(10.0 ** (-q / 10.0))))


def read_length(read: pysam.AlignedSegment) -> int:
    n = read.query_length
    if not n:
        n = read.infer_read_length() or 0
    return int(n)


def read_label(read: pysam.AlignedSegment, bam_path: str, label_by: str = "auto") -> str:
    """Sample label for per-sample thresholds / plots.

    label_by:
      auto  - BC tag barcode (e.g. 'barcode09'), else barcodeNN in the file
              name, else the file stem
      file  - file stem
      bc    - BC tag (full value) or 'nobc'
      rg    - RG tag or 'norg'
      rg+bc - '<RG>|<BC>' (safe for merged multi-run BAMs)
      none  - 'all'
    """
    if label_by == "none":
        return "all"
    stem = os.path.basename(bam_path)
    stem = stem[:-4] if stem.endswith(".bam") else stem
    if label_by == "file":
        return stem
    bc = read.get_tag("BC") if read.has_tag("BC") else None
    rg = read.get_tag("RG") if read.has_tag("RG") else None
    if label_by == "bc":
        return str(bc) if bc else "nobc"
    if label_by == "rg":
        return str(rg) if rg else "norg"
    if label_by == "rg+bc":
        return f"{rg}|{bc}"
    # auto
    if bc:
        m = _BARCODE_RE.search(str(bc))
        return m.group(1).lower() if m else str(bc)
    m = _BARCODE_RE.search(stem)
    return m.group(1).lower() if m else stem


def mm_mode(read: pysam.AlignedSegment) -> str:
    """MM mode flag of the m6A group: '.', '?' or '' (unspecified = '.')."""
    for t in ("MM", "Mm"):
        if read.has_tag(t):
            tag = str(read.get_tag(t))
            for group in tag.split(";"):
                head = group.split(",")[0]
                if head.startswith(("A+a", "A-a", "A+21891", "A-21891")):
                    return head[-1] if head[-1] in ".?" else "."
            return "."
    return "."


def m6a_ml(read: pysam.AlignedSegment) -> Tuple[np.ndarray, str]:
    """Dense ML array over stored-query positions (-1 = not listed in MM).

    Returns (ml int16[qlen], mm_mode).
    """
    seq = read.query_sequence
    qlen = len(seq) if seq else 0
    ml = np.full(qlen, -1, dtype=np.int16)
    mb = read.modified_bases
    if mb:
        for k, v in mb.items():
            if k[0] == "A" and k[2] in M6A_CODES and len(v):
                arr = np.asarray(v, dtype=np.int64)
                ml[arr[:, 0]] = arr[:, 1]
    return ml, mm_mode(read)


def target_mask(read: pysam.AlignedSegment) -> np.ndarray:
    """Boolean over stored SEQ: True where the base is the m6A target
    (A on forward reads, T on reverse reads = template A)."""
    seq = read.query_sequence or ""
    tb = b"T" if read.is_reverse else b"A"
    return np.frombuffer(seq.encode(), dtype="S1") == tb


def aligned_query_ref(read: pysam.AlignedSegment) -> Tuple[np.ndarray, np.ndarray]:
    """Vectorised get_aligned_pairs(matches_only=True): (qpos[], rpos[])."""
    qs, rs = [], []
    q = 0
    r = read.reference_start
    for op, n in read.cigartuples or ():
        if op in (0, 7, 8):  # M, =, X
            qs.append(np.arange(q, q + n, dtype=np.int64))
            rs.append(np.arange(r, r + n, dtype=np.int64))
            q += n
            r += n
        elif op in (1, 4):   # I, S
            q += n
        elif op in (2, 3):   # D, N
            r += n
        # H, P: consume nothing
    if not qs:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(qs), np.concatenate(rs)


# --------------------------------------------------------------------------
# thresholds
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class CallThreshold:
    ml_meth: int
    ml_canon: int

    @classmethod
    def fixed(cls, ml: int) -> "CallThreshold":
        """Binary call: ML >= ml is methylated, everything else canonical."""
        return cls(int(ml), int(ml) - 1)

    @classmethod
    def from_conf(cls, thr: float) -> "CallThreshold":
        """modkit-style confidence threshold -> integer ML cut-offs."""
        ml_meth = int(round(thr * 256 - 0.5))
        return cls(ml_meth, 255 - ml_meth)

    @property
    def conf(self) -> float:
        return (self.ml_meth + 0.5) / 256.0


def classify(ml: np.ndarray, mode: str, thr: CallThreshold) -> np.ndarray:
    """ML array -> call array (1 meth, 0 canonical, -1 no-call)."""
    calls = np.full(ml.shape, -1, dtype=np.int8)
    listed = ml >= 0
    calls[listed & (ml >= thr.ml_meth)] = 1
    calls[listed & (ml <= thr.ml_canon)] = 0
    if mode != "?":
        calls[~listed] = 0
    return calls


# --------------------------------------------------------------------------
# fibers
# --------------------------------------------------------------------------
@dataclass
class Fiber:
    read_id: str
    chrom: str
    start: int            # reference_start (0-based)
    end: int              # reference_end (exclusive)
    strand: str           # '+' / '-'
    label: str
    pos: np.ndarray       # int64 ref positions of aligned target A's
    call: np.ndarray      # int8 call state at pos
    mapq: int = 0
    extra: dict = field(default_factory=dict)

    @property
    def meth_pos(self) -> np.ndarray:
        return self.pos[self.call == 1]

    @property
    def n_meth(self) -> int:
        return int((self.call == 1).sum())

    @property
    def n_called(self) -> int:
        return int((self.call >= 0).sum())

    @property
    def rate(self) -> float:
        n = self.n_called
        return self.n_meth / n if n else float("nan")


def read_to_fiber(read: pysam.AlignedSegment, thr: CallThreshold, label: str,
                  require_ref_match: bool = False,
                  ref_seq: Optional[str] = None, ref_offset: int = 0) -> Fiber:
    """Make a Fiber (calls at aligned target-A reference positions)."""
    ml, mode = m6a_ml(read)
    calls_q = classify(ml, mode, thr)
    tmask = target_mask(read)
    qpos, rpos = aligned_query_ref(read)
    keep = tmask[qpos]
    if require_ref_match and ref_seq is not None:
        tb = "T" if read.is_reverse else "A"
        idx = rpos - ref_offset
        ok = (idx >= 0) & (idx < len(ref_seq))
        refb = np.array([ref_seq[i] if o else "N" for i, o in zip(idx, ok)])
        keep &= np.char.upper(refb) == tb
    qpos, rpos = qpos[keep], rpos[keep]
    return Fiber(read.query_name, read.reference_name, read.reference_start,
                 read.reference_end, "-" if read.is_reverse else "+", label,
                 rpos, calls_q[qpos], read.mapping_quality)


def read_space_counts(read: pysam.AlignedSegment, thr: CallThreshold) -> Tuple[int, int, int]:
    """Counts over ALL target bases in the stored read (incl. soft clips):
    (n_meth, n_canonical, n_target)."""
    ml, mode = m6a_ml(read)
    c = classify(ml, mode, thr)[target_mask(read)]
    return int((c == 1).sum()), int((c == 0).sum()), int(c.size)


class Thresholds:
    """Per-label call thresholds with a default.

    Build with Thresholds.fixed(192), Thresholds.from_tsv(path) or
    fibernhat.io.thresholds.estimate(...).
    """

    def __init__(self, by_label: Optional[Dict[str, CallThreshold]] = None,
                 default: Optional[CallThreshold] = None, source: str = ""):
        self.by_label = dict(by_label or {})
        self.default = default
        self.source = source

    @classmethod
    def fixed(cls, ml: int) -> "Thresholds":
        return cls({}, CallThreshold.fixed(ml), source=f"fixed ML>={ml}")

    def get(self, label: str) -> CallThreshold:
        if label in self.by_label:
            return self.by_label[label]
        if self.default is not None:
            return self.default
        raise KeyError(f"no m6A threshold for label {label!r}; known: "
                       f"{sorted(self.by_label)}. Run `fibernhat thresholds` "
                       f"or pass --threshold-ml.")

    def __repr__(self):
        return f"Thresholds({self.source}; {len(self.by_label)} labels; default={self.default})"


def open_bam(path: str) -> pysam.AlignmentFile:
    return pysam.AlignmentFile(path, "rb", check_sq=False)


def iter_reads(bam_paths: Sequence[str], region: Optional[Tuple[str, int, int]] = None,
               filt: Optional[ReadFilter] = None, label_by: str = "auto",
               ) -> Iterator[Tuple[pysam.AlignedSegment, str, str]]:
    """Yield (read, label, bam_path) for reads passing `filt` (None = no filter
    except primary mapped)."""
    filt = filt or ReadFilter.none()
    for p in bam_paths:
        with open_bam(p) as bam:
            if region is not None:
                it = bam.fetch(region[0], region[1], region[2])
            else:
                it = bam.fetch(until_eof=True)
            for read in it:
                if filt.reason(read) is not None:
                    continue
                yield read, read_label(read, p, label_by), p


def fetch_fibers(bam_paths: Sequence[str], region: Tuple[str, int, int],
                 thresholds: Thresholds, filt: Optional[ReadFilter] = None,
                 label_by: str = "auto", require_span: bool = False,
                 max_reads: Optional[int] = None, seed: int = 0) -> List[Fiber]:
    """All fibers overlapping `region` (calls restricted to the region)."""
    chrom, s, e = region
    fibers = []
    for read, label, _ in iter_reads(bam_paths, region, filt, label_by):
        if require_span and not (read.reference_start <= s and read.reference_end >= e):
            continue
        f = read_to_fiber(read, thresholds.get(label), label)
        m = (f.pos >= s) & (f.pos < e)
        f.pos, f.call = f.pos[m], f.call[m]
        fibers.append(f)
    if max_reads is not None and len(fibers) > max_reads:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(len(fibers), max_reads, replace=False))
        fibers = [fibers[i] for i in idx]
    return fibers
