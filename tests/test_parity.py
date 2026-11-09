"""pysam engine == modkit engine (modkit 0.4.4). Required by the repo rules."""
import numpy as np

from fibernhat import pileup as P
from fibernhat import quintiles as Q
from fibernhat.io import thresholds as T
from fibernhat.io.bam import CallThreshold, ReadFilter, Thresholds


def _th(df):
    return Thresholds({r.label: CallThreshold(int(r.ml_meth), int(r.ml_canon)) for r in df.itertuples()})


def test_threshold_parity(bams, modkit):
    for filt in (ReadFilter(), ReadFilter.none()):
        a = T.estimate(bams, filt, "file", max_reads_per_label=None)
        b = T.estimate(bams, filt, "file", max_reads_per_label=None, engine="modkit", modkit=modkit)
        assert a.ml_meth.tolist() == b.ml_meth.tolist()
        assert np.allclose(a.thr, b.thr)


def test_pileup_parity_per_label(bams, modkit):
    th = _th(T.estimate(bams, ReadFilter(), "file", max_reads_per_label=None))
    d = P.compare(P.pileup(bams, th, ReadFilter(), "file"),
                  P.pileup(bams, th, ReadFilter(), "file", engine="modkit", modkit=modkit))
    assert len(d) == 0, d.head()


def test_pileup_parity_fixed128(bams, modkit):
    th = Thresholds.fixed(128)
    d = P.compare(P.pileup(bams, th, ReadFilter(), "file"),
                  P.pileup(bams, th, ReadFilter(), "file", engine="modkit", modkit=modkit))
    assert len(d) == 0, d.head()


def test_quintile_parity(bams, modkit):
    th = _th(T.estimate(bams, ReadFilter(), "auto", max_reads_per_label=None))
    _, p1, s1 = Q.run(bams, th, ReadFilter(), ("dyad", "tss"), progress=False)
    _, p2, s2 = Q.run(bams, th, ReadFilter(), ("dyad", "tss"), engine="modkit", modkit=modkit,
                      progress=False)
    assert (p1.n_mod.values == p2.n_mod.values).all()
    assert (p1.n_valid.values == p2.n_valid.values).all()
    assert (s1.genome_n_valid.values == s2.genome_n_valid.values).all()
