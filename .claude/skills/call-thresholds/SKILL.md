---
name: call-thresholds
description: Compute per-sample (per-barcode) m6A call thresholds from BAM — the 10th-percentile-confidence cutoff (modkit sample-probs semantics) used by every fibernhat plot. Use when asked about m6A thresholds, ML cutoffs, call confidence, or to make plots reproducible across runs.
---
# Per-sample m6A call thresholds

## Rule
For each target base (A on + reads, T on − reads, soft clips included):
- p = (ML + 0.5)/256 and conf = max(p, 1 − p). A's missing from MM count as canonical with conf 1.0 when the MM mode is `.`.
- **thr** = the 10th percentile of conf. Then **ml_meth** = round(thr·256 − 0.5) and **ml_canon** = 255 − ml_meth.
- Call: methylated if ML ≥ ml_meth; canonical if ML ≤ ml_canon or the A is unlisted; anything in between is a **no-call** and is ignored.

## Command
```bash
.venv/bin/fibernhat thresholds a.bam b.bam -o outputs/m6a_thresholds.tsv            # 50k reads/sample
.venv/bin/fibernhat thresholds a.bam b.bam -o thr.tsv --all-reads                    # like the lab TSV
.venv/bin/fibernhat thresholds a.bam -o thr.tsv --engine modkit --all-reads          # same numbers, via modkit
```
Pass the result to any plotting command with `--thresholds thr.tsv`. The lab TSV `ref/thresholds/call_thresholds_by_barcode_filtnordna.tsv` also works: it matches reads by their RG.

## Notes
- The default read filters are applied before estimating, the same as the lab TSV, which was computed on `filtered.bam`.
- With the modkit engine, each sample is written to a filtered BAM and `modkit sample-probs --no-sampling` is run on it. pysam and modkit give **identical** thresholds; this is checked in `tests/test_parity.py`.
- Typical values: Hia5 samples are about 0.78–0.83 (ML 199–212). no-MTase controls are about 0.90.
- If any plotting command runs without `--thresholds` or `--threshold-ml`, it estimates thresholds on the fly from 20k reads per sample and prints them to stderr.

## Provenance
roboNhat `analysis/singlefiber/ml_histograms.py` (`modkit_conf`, `percentile_thr`) and `read_calls.py` (call rule). The lab TSV was made by `filtered_refit_nordna_bcthr/scripts/job1_split_sample.sh`.
