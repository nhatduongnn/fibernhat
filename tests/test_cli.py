"""Smoke tests: every command runs on the test BAMs and writes its outputs."""
import os

import pandas as pd

from fibernhat.cli import main


def test_qc(tmp_path, bams):
    main(["qc", *bams, "-o", str(tmp_path), "--formats", "png"])
    s = pd.read_csv(tmp_path / "summary.tsv", sep="\t")
    assert set(s.label) == {"barcode09", "barcode10"}
    assert (s.n50 > 0).all()
    assert (tmp_path / "summary_read_length.png").exists()


def test_fiberplot(tmp_path, bams):
    main(["fiberplot", "--bam", *bams, "--region", "chrII:277500-280000", "--threshold-ml", "192",
          "--no-mnase", "-o", str(tmp_path), "--formats", "png"])
    assert (tmp_path / "chrII_277500-280000.png").exists()
    t = pd.read_csv(tmp_path / "chrII_277500-280000_fibers.tsv", sep="\t")
    assert len(t) > 50


def test_aggregate(tmp_path, bams):
    main(["aggregate", "--bam", *bams, "--sites", "tss", "--chrom", "chrI", "--half", "300",
          "--no-mnase", "--threshold-ml", "192", "-o", str(tmp_path), "--formats", "png"])
    p = pd.read_csv(tmp_path / "tss_profile.tsv", sep="\t")
    assert p.n_valid.sum() > 0


def test_quintiles(tmp_path, bams):
    main(["quintiles", "--bam", *bams, "--sites", "dyad,tss", "-o", str(tmp_path), "--formats", "png"])
    s = pd.read_csv(tmp_path / "quintile_summary.tsv", sep="\t")
    # quintile genome m6A rate increases monotonically Q1 -> Q5 within each sample
    for _, g in s.groupby("label"):
        assert g.sort_values("quintile").genome_m6a_rate.is_monotonic_increasing


def test_fiberplot_gene(tmp_path, bams):
    main(["fiberplot", "--bam", *bams, "--gene", "TEA1", "--flank", "500", "--threshold-ml", "192",
          "--no-mnase", "-o", str(tmp_path), "--formats", "png"])
    assert (tmp_path / "TEA1.png").exists()


def test_panel_order(tmp_path, bams):
    import pytest
    from fibernhat.plot.fibers import parse_panels
    assert parse_panels("genes,mnase,fibers,pileup") == ["genes", "mnase", "fibers", "pileup"]
    with pytest.raises(ValueError):
        parse_panels("genes,heatmap")
    main(["fiberplot", "--bam", *bams, "--region", "chrII:277500-280000", "--threshold-ml", "192",
          "--no-mnase", "--panels", "pileup,fibers,genes", "--pool", "--name", "reordered",
          "-o", str(tmp_path), "--formats", "png"])
    assert (tmp_path / "reordered.png").exists()
    main(["aggregate", "--bam", *bams, "--sites", "tss", "--chrom", "chrI", "--half", "300",
          "--no-mnase", "--threshold-ml", "192", "--panels", "ratio,fibers", "--name", "agg2",
          "-o", str(tmp_path), "--formats", "png"])
    assert (tmp_path / "agg2.png").exists()
