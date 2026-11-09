# Capabilities

One row per skill. Status values: **stable** (tested, parity-checked where it applies), **beta** (works, may still change), **planned** (code exists elsewhere in the lab and has not been ported yet).

| skill | status | command | input | output | ported from |
|---|---|---|---|---|---|
| [`read-qc`](.claude/skills/read-qc/SKILL.md) | stable | `fibernhat qc` | any BAM(s) | length/Q histograms, violins, N50/coverage table, filter cascade | `process_nanopore_sequencing/Nov_20_25/no_trim/plot_read_data.py` |
| [`call-thresholds`](.claude/skills/call-thresholds/SKILL.md) | stable, parity | `fibernhat thresholds` | BAM(s) with MM/ML | per-sample 10th-pct-confidence ML cutoffs (TSV) | roboNhat `singlefiber/ml_histograms.py`, modkit `sample-probs` |
| [`raw-fiberplot`](.claude/skills/raw-fiberplot/SKILL.md) | stable | `fibernhat fiberplot` | BAM(s) + region | gene track, MNase V-plot, per-read dots, m6A pileup (panel order configurable) | `hmm…/scripts/fig2_panels.py`, `fig1_panels.py`, Strahl_lab TEA1, DMS_seq annotation |
| [`aggregate-fiberplot`](.claude/skills/aggregate-fiberplot/SKILL.md) | beta | `fibernhat aggregate` | BAM(s) + site set (dyad/tss/abf1/reb1/tf/file) | piled-up fibers, m6A profile, MNase V-plot, npz cache | `Fiber_seq_plots.process_aggregate_fiber_plot`, `aggregate_fiberplot_open_sites.py` |
| [`quintile-metaprofile`](.claude/skills/quintile-metaprofile/SKILL.md) | stable, parity | `fibernhat quintiles` | BAM(s) | per-read m6A-rate quintile profiles at dyad/TSS/Abf1/Reb1, footprint scores | `investigate_un_crosslinked_runs/per_read_heterogeneity.py`, `footprint_metaprofiles.py` |
| [`basecall-align-dorado`](.claude/skills/basecall-align-dorado/SKILL.md) | beta | `sbatch workflows/slurm/basecall_align_demux.sbatch` | pod5 | aligned demuxed BAMs + QC + thresholds | `process_nanopore_sequencing/DM1746/script_sup_model.sh` |
| pileup (utility) | stable, parity | `fibernhat pileup`, `fibernhat parity` | BAM(s) | per-position n_mod/n_valid (bedMethyl-like); pysam==modkit check | — |

## Planned (not yet ported)

| skill | where the code lives now | notes |
|---|---|---|
| `decode-hsmm` | `projects/Fiber_seq/hmm_for_clustering_claude_code` (`scripts/decode_chrom_save.py`, `src/models/kmer_hsmm.py`, model `output/erv46_4st_10iter_flip/model.npz`, 48 KB) | 4-state k-mer HSMM (open/linker/TF/nuc); posterior P(accessible) |
| `cluster-fibers` | same repo (`extract_fiber_features.py`, `cluster_fibers.py`), `collab/Strahl_lab/scripts/build_tea1_plots.py`, `relabel_clusters_by_nuc.py` | UMAP+HDBSCAN or Ward on posterior; clustered raw/posterior fiberplots |
| `persite-quintiles` | `hmm…/scripts/fig5_build_persite_quintile_cache.py`, `fig5_render_persite_quintile.py` | fibers ranked per site (fragile nucleosomes) |
| `robocop-permol` | `programs/roboNhat_w_new_changes_claude/analysis/singlefiber/` | per-read RoboCOP decode, needs C library and a 94 MB HMMconfig |
| `ml-histograms` | roboNhat `singlefiber/ml_histograms.py` | ML distributions by region (linker / core / TF) |
| `tf-pileup` | roboNhat `pileup_tf_windows.py`, `plot_tf_pileup.py` | per-TF Watson/Crick profiles |

To add a skill, see "Adding a new skill" in `AGENTS.md`.
