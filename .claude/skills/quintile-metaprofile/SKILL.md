---
name: quintile-metaprofile
description: Split reads into quintiles by their own per-read m6A rate and plot m6A metaprofiles (one line per quintile) around nucleosome dyads, TSSs and Abf1/Reb1 sites; also footprint scores per quintile. Use for "quintile plot", "per-read methylation heterogeneity", "signal around nucs/TF/TSS by m6A quintile".
---
# Per-read m6A-rate quintile metaprofiles

## What it does
1. **Pass 1:** for each read that passes the filters and is ≥ `--min-len` (500) bp on a nuclear chromosome, compute rate = methylated / called A over the whole read. Then, separately for each sample, `qcut` the rates into 5 groups: Q1 = least methylated, Q5 = most.
2. **Pass 2:** for each quintile, build a per-position m6A pileup (both strands summed) and sum it over oriented site windows. The windows are dyad ±1000, TSS ±1500, Abf1/Reb1 ±300. Each quintile is also given its genome-wide m6A rate.
3. **Plot:** one row per sample and one column per site class. Each line is m6A divided by that quintile's genome mean, smoothed over 21 bp, coloured with viridis.

## Command
```bash
.venv/bin/fibernhat quintiles --bam a.bam b.bam -o outputs/quintiles                 # dyad,tss,abf1,reb1
.venv/bin/fibernhat quintiles --bam a.bam --sites dyad,tss,chereji --half dyad=800
# reproduce the original Step-6 analysis exactly:
.venv/bin/fibernhat quintiles --bam a.bam --threshold-ml 128 --min-mapq 10 --min-qs none --keep-rdna
# same result through modkit pileup (needs modkit 0.4.4):
.venv/bin/fibernhat quintiles --bam a.bam --engine modkit
```

## Outputs
- `quintile_profiles.png/pdf`
- `per_read_rates.tsv.gz`: rate and quintile for each read
- `metaprofiles_raw.tsv.gz`: label, quintile, site_class, pos, n_mod, n_valid, m6a_rate, n_sites
- `quintile_summary.tsv`: n_reads, rate range, genome m6A rate per quintile
- `footprint_scores.tsv`: `dyad_linker_over_dyad`, `dyad_nrl_power`, `tss_ndr_over_body`, `<tf>_flank_over_center`

## Notes
- These quintiles group **reads by their own m6A rate**. They are not expression quintiles, and they are not the per-site fibre ranking of figure 5.
- The pysam and modkit engines give identical profiles (`tests/test_parity.py`).
- The original analysis subsampled every sample to the same read count. Do that beforehand with `samtools view -s` if you need it.
- A full run per barcode takes about 1–3 passes over the BAM. Use `workflows/slurm/fibernhat_job.sbatch` for big inputs.

## Provenance
`Fiber_seq/investigate_un_crosslinked_runs/per_read_heterogeneity.py` (`per_read_rates`, `stage_rates`, `plot_quintiles`), `footprint_metaprofiles.py` (`load_sites`, `window_sum`, `score_sample`), `plot_quintile_pair.py`, and `scripts/06{a,b,c}_*.sbatch`. The figure it reproduces is `results/step6/quintile_profiles_Mar20_25_bc01_vs_Mar20_25_bc03.png`.
