"""SGD ORF / Park 2014 TSS-PAS gene annotation track (vendored DMS_seq plotter)."""
from __future__ import annotations

import logging
from functools import lru_cache

from ..io import refs

logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

_ROMAN_TO_INT = {r: i + 1 for i, r in enumerate(refs.ROMAN)}


def chrom_str_to_int(chrom: str):
    """'chrXVI' -> 16; None for chromosomes without SGD annotation (chrM ...)."""
    name = chrom[3:] if chrom.startswith("chr") else chrom
    return _ROMAN_TO_INT.get(name)


@lru_cache(maxsize=4)
def make_orf_annotation_plotter(sgd_gff: str = str(refs.SGD_GFF),
                                tss_gff: str = str(refs.TSS_GFF),
                                pas_gff: str = str(refs.PAS_GFF)):
    from .reference_data import read_sgd_orfs, read_sgd_orf_introns, read_park_TSS_PAS
    from .orf_annotation_plotter import ORFAnnotationPlotter
    orfs = read_sgd_orfs(sgd_gff).join(read_park_TSS_PAS(tss_gff, pas_gff))
    introns = read_sgd_orf_introns(sgd_gff)
    return ORFAnnotationPlotter(orfs, introns)


def find_gene(name: str):
    """Look up a gene by common or systematic name (case-insensitive).

    Returns dict(name, orf, chrom, start, stop (1-based), strand, tss, tss_source):
    tss = Park 2014 TSS if present, else the start codon (tss_source says which).
    """
    import math
    p = make_orf_annotation_plotter()
    orfs = p.orfs if hasattr(p, "orfs") else None
    if orfs is None:
        from .reference_data import read_sgd_orfs, read_park_TSS_PAS
        orfs = read_sgd_orfs(str(refs.SGD_GFF)).join(read_park_TSS_PAS(str(refs.TSS_GFF), str(refs.PAS_GFF)))
    key = name.upper()
    hit = orfs[(orfs.name.astype(str).str.upper() == key) | (orfs.index.astype(str).str.upper() == key)]
    if hit.empty:
        raise ValueError(f"gene {name!r} not found in SGD annotation")
    r = hit.iloc[0]
    chrom = f"chr{refs.ROMAN[int(r.chr) - 1]}"
    start, stop = int(min(r.start, r.stop)), int(max(r.start, r.stop))
    atg = start if r.strand == "+" else stop
    tss = r.get("TSS")
    has_tss = tss is not None and not (isinstance(tss, float) and math.isnan(tss))
    return dict(name=str(r["name"]), orf=str(hit.index[0]), chrom=chrom, start=start, stop=stop,
                strand=r.strand, tss=int(tss) if has_tss else atg,
                tss_source="Park 2014 TSS" if has_tss else "start codon (no Park TSS)")
