import numpy as np
import pandas as pd

from fibernhat.io import refs
from fibernhat.io import thresholds as T
from fibernhat.io.bam import (CallThreshold, ReadFilter, Thresholds, classify, fetch_fibers,
                              open_bam, parse_region, read_label, read_to_fiber, target_mask)
from fibernhat.qc.read_stats import n50


def test_parse_region():
    assert parse_region("chrI:60,500-63,500") == ("chrI", 60499, 63500)


def test_classify_modes():
    ml = np.array([-1, 0, 54, 55, 100, 199, 200, 255], np.int16)
    thr = CallThreshold(200, 55)
    assert classify(ml, ".", thr).tolist() == [0, 0, 0, 0, -1, -1, 1, 1]
    assert classify(ml, "?", thr).tolist() == [-1, 0, 0, 0, -1, -1, 1, 1]
    fx = CallThreshold.fixed(128)
    assert classify(np.array([127, 128, -1], np.int16), ".", fx).tolist() == [0, 1, 0]


def test_conf_roundtrip():
    for k in range(128, 256):
        thr = (k + 0.5) / 256
        ct = CallThreshold.from_conf(thr)
        assert ct.ml_meth == k and ct.ml_canon == 255 - k


def test_percentile_thr_simple():
    h = np.zeros(257, np.int64)
    h[200] = 10          # conf (200.5/256)
    h[256] = 90          # unlisted, conf 1.0
    assert T.percentile_thr(h, 0.10) == (200 + 0.5) / 256


def test_n50():
    assert n50([1, 2, 3, 4, 10]) == 10
    assert n50([5, 5, 5, 5]) == 5


def test_labels(bams):
    with open_bam(bams[0]) as b:
        r = next(b.fetch("chrI"))
    assert read_label(r, bams[0], "auto") == "barcode09"
    assert read_label(r, bams[0], "file") == "nov20_barcode09_subset"


def test_reverse_reads_scored_at_ref_T(bams):
    fa = __import__("pysam").FastaFile(str(refs.FASTA))
    n_rev = 0
    with open_bam(bams[0]) as b:
        for r in b.fetch("chrI", 60000, 64000):
            if not ReadFilter().passes(r):
                continue
            f = read_to_fiber(r, CallThreshold.fixed(128), "x")
            if f.pos.size == 0:
                continue
            ref = np.array([fa.fetch("chrI", p, p + 1).upper() for p in f.pos[:200]])
            want = "T" if r.is_reverse else "A"
            assert (ref == want).mean() > 0.9   # aligned target bases sit on ref A/T (mismatches allowed)
            n_rev += r.is_reverse
    assert n_rev > 0


def test_fetch_fibers_region(bams):
    th = Thresholds.fixed(192)
    fb = fetch_fibers(bams, parse_region("chrII:277500-280000"), th)
    assert len(fb) > 50
    for f in fb:
        assert ((f.pos >= 277499) & (f.pos < 280000)).all()


def test_site_loaders():
    assert len(refs.load_sites("dyad")) == 2000
    t = refs.load_sites("tss")
    assert len(t) > 5000 and set(t.sign) == {-1, 1}
    a = refs.load_sites("abf1")
    assert a.chrom.str.startswith("chr").all() and "chr9" not in set(a.chrom)
    assert refs.nuclear_genome_size() == 12071326


def test_thresholds_tsv_roundtrip(tmp_path, bams):
    df = T.estimate(bams, ReadFilter(), "auto", max_reads_per_label=300)
    p = tmp_path / "t.tsv"
    T.write_tsv(df, str(p))
    th = T.read_tsv(str(p))
    for r in df.itertuples():
        assert th.get(r.label).ml_meth == r.ml_meth


def test_lab_tsv_loads():
    th = T.read_tsv(str(refs.REF_DIR / "thresholds" / "call_thresholds_by_barcode_filtnordna.tsv"))
    assert th.get("Mar20_bc01").ml_meth == 210


def test_find_gene():
    from fibernhat.annotation import find_gene
    g = find_gene("tea1")
    assert g["chrom"] == "chrXV" and g["tss"] == 954213 and g["tss_source"].startswith("Park")
    assert find_gene("YBR020W")["name"] == "GAL1"
