# Reference files

These are all small, so they're committed to git. Each file below lists where it came from on the lab cluster (Duke CS `/usr/xtmp/nd141`).

## `sacCer3/`

| file | what | source |
|---|---|---|
| `custom_sacCer3_lambda.fasta` (+`.fai`) | sacCer3 (UCSC names `chrI`…`chrXVI`, `chrM`) **plus `chrLambda` (3,560 bp)**. All 17 yeast contigs are identical to `yeast_UCSC-sacCer3_SGD-R64.fasta`. | `projects/sequence_alignment/data/ref-data/yeast/` |
| `sacCer3_lambda.chrom.sizes` | contig lengths that match the FASTA | `.../yeast/yeast_genome_w_Lambda.genome` |
| `saccharomyces_cerevisiae_R64-1-1_20110208_no_fasta.gff` | SGD R64-1-1 ORFs/introns (gene track) | `projects/DMS_seq/src/script/` |
| `Park_2014_TSS_V64.gff`, `Park_2014_PAS_V64.gff` | Park et al. 2014 TSS / PAS (gene track arrows, TSS site set) | `projects/DMS_seq/src/script/` |

We picked the lambda FASTA because it is what the lab's Nov_20_25 and DM1746 runs were aligned to: their BAM `@SQ` lines (including `M5` md5s) match it exactly. BAMs aligned to plain sacCer3 also work, because contigs are matched by name.

**Do not use `yeast_ucsc_saccer3_yam095.fasta`** (strain yAM095). Its chrIII has a ~97 kb segment between MAT and HMR replaced (+1,265 bp, including vector sequence), so its coordinates don't match standard sacCer3.

## `sites/`

| file | n | used as | source |
|---|---|---|---|
| `Brogaard_nuc_positions.sacCer3.top2000.tsv` | 2000 | `--sites dyad` (1-based) | `projects/replicate_prob_dyad_plot/data/ref-data/` |
| `Chereji_2018_+1_-1_nucs.bed` | 9271 | `--sites chereji` | roboNhat `analysis/inputs/` |
| `abf1_fimo.bed`, `reb1_fimo.bed` | 569 / 698 with motif | `--sites abf1/reb1` (motif midpoint, motif strand) | `projects/data/lab_shared_dataset/rossi_strand/output/` |
| `rossi_peak_w_strand_conformed_to_PWM_all_TFs_peakVal.bed` | 6000 | `--sites tf` / `tf:<TF>` | roboNhat `analysis/inputs/` |
| `MacIsaac_sacCer3_liftOver_Abf1_Reb1.bed` | 307 | any `--sites <path>` | roboNhat `analysis/inputs/` |
| `park_tss_sites.tsv`, `abf1_rossi_sites.tsv` | 5203 / 569 | paper figure 1 site tables | `projects/Fiber_seq/paper_draft/figure_1/` |

**Note on the dyad file:** the original quintile analysis (`investigate_un_crosslinked_runs/footprint_metaprofiles.py`) used `Brogaard_nuc_positions.top2000.tsv`. That file is the same list but some positions are shifted by about 2 bp. We ship the `sacCer3`-named version.

## `thresholds/`

`call_thresholds_by_barcode_filtnordna.tsv` holds the lab's per-barcode m6A cutoffs for the 7 barcodes in the merged training BAM. They were made with modkit 0.2.6 `sample-probs` and the 10th percentile. `fibernhat --thresholds` reads this file directly, matching reads by their RG or barcode name.

## md5

```
305e65e553371d0660d97b204b6b9b34 sacCer3/custom_sacCer3_lambda.fasta
37274fb87e4f6ed384c9e3f8294d9bcc sacCer3/saccharomyces_cerevisiae_R64-1-1_20110208_no_fasta.gff
499fe43c006081e4beea7734a8fffb29 sacCer3/Park_2014_TSS_V64.gff
fb2c337e25220fd2a9c85ced06fcf568 sacCer3/Park_2014_PAS_V64.gff
8602837229f2cc24c71b37fac41faab2 sites/Brogaard_nuc_positions.sacCer3.top2000.tsv
```
