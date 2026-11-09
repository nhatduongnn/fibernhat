---
name: aggregate-fiberplot
description: Aggregate (piled-up) fiber-seq plots around many sites — nucleosome dyads, TSSs, TF motifs (Abf1/Reb1/any Rossi TF) or a custom BED/TSV — with per-fiber dot rows, m6A/called-A profile and optional pooled MNase V-plot. Use for "aggregate fiberplot", "metaplot around nucleosomes/TSS/TF", "footprint at Abf1".
---
# Aggregate fiberplot around sites

## Command
```bash
.venv/bin/fibernhat aggregate --bam a.bam b.bam --sites dyad -o outputs/agg
.venv/bin/fibernhat aggregate --bam a.bam --sites tss --half 1000 --sort meth_in_win
.venv/bin/fibernhat aggregate --bam a.bam --sites tf:Reb1 --half 300
.venv/bin/fibernhat aggregate --bam a.bam --sites my_sites.tsv --half 500   # chrom,pos[,strand] (1-based) or BED
```
Built-in site sets (see `ref/README.md`):
- `dyad` (Brogaard top 2000)
- `chereji` (Chereji ±1 nucleosomes)
- `tss` (Park 2014, oriented)
- `abf1`, `reb1` (Rossi ChIP-exo summits with a FIMO motif, oriented)
- `tf` / `tf:<Name>` (Rossi PWM-conformed peaks)

## Outputs
- `<name>.png/pdf`, top to bottom (default `--panels mnase,fibers,counts,ratio`): MNase V-plot (if an MNase BAM is configured), one dot panel per sample, "Methylation pileup" counts, then "Methylated As / Possible As". Counts and ratio are unsmoothed, as in `process_aggregate_fiber_plot`; `--smooth N` adds an N bp rolling mean.
- `<name>_profile.tsv`: per sample and position, n_mod, n_valid and m6a_frac.
- `<name>_fibers.npz`: every (fiber, site) row with oriented methylated/called positions. Reload with `fibernhat.aggregate.SiteFibers.load`.

## Options
- `--panels`: any subset and order of `mnase,fibers,counts,ratio`.
- `--sort`:
  - `meth_in_win` (default, as in `aggregate_fiberplot_open_sites.py`): m6A fraction within ±100 bp, most methylated on top
  - `shuffle`
  - `meth_frac`: whole window
  - `site`
- `--max-rows` caps how many rows are drawn per sample. The profile always uses every fiber.
- `--max-per-site`, `--max-sites`, `--chrom chrI,chrII` let you take a quick subset.
- By default a read must span the whole window. `--no-require-span` includes reads that only cover part of it.
- `--no-orient` stops minus-strand sites from being flipped.

## Notes
- Strand-aware: windows are flipped for − sites so the site's 5'→3' direction points right.
- MNase uses `mnase_bam` from the site config (fragments ≤ 250 bp, up to 10k points).
- Genome-wide site sets on a full BAM take minutes. For many BAMs, use `workflows/slurm/fibernhat_job.sbatch`.

## Provenance
Pure-numpy reimplementation of `Fiber_seq/src/script/Fiber_seq_plots.py::process_aggregate_fiber_plot` and `hmm_for_clustering_claude_code/scripts/aggregate_fiberplot_open_sites.py` (`WindowedDecodeMixin`, `stamp_meth_in_win`), with no pyranges or parquet. The look follows `fig1_panels.draw_fiber_heatmap` and `draw_mnase_vplot`.
