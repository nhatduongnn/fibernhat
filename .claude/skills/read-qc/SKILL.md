---
name: read-qc
description: Plot read length and read quality QC (N50, histograms, violins, filter cascade) for nanopore fiber-seq BAMs. Use when someone wants QC, read length, quality, N50, coverage or a "how did the run go" summary from BAM files (aligned or unaligned, any dorado version, any cluster).
---
# Read QC from BAM

**When to use:** QC of a sequencing run or barcode. You only need the BAM: it can be aligned or unaligned, one per barcode or merged. Dorado output and the lab's directory layout are not required.

## Command
```bash
.venv/bin/fibernhat qc path/to/*.bam -o outputs/qc_<run>
# quick look: --max-reads 50000 ; plot only reads passing filters: --pass-only
# merged multi-run BAM: --label-by rg+bc ; one sample per file: --label-by file
```

## Outputs (in `-o`)
| file | content |
|---|---|
| `summary.tsv` | per sample: n_reads, total_bases, **N50**, median/mean length, median/mean Q, frac ≥ Q10, frac mapped, frac with MM tag, coverage, reads and coverage passing the filters |
| `read_stats.tsv.gz` | one row per primary read: length, mean Q, MAPQ, ref span, `fail_reason` |
| `filter_cascade.tsv/.png` | reads removed at each filter step (unmapped, no MM, chrM/lambda, rDNA, MAPQ<20, span<500, qs<10) |
| `<label>_hist.png/.svg` | length histogram weighted by bases (10 bp bins, clipped at `--clip`), KDE, N50 and median lines |
| `summary_read_length.*`, `summary_quality.*` | violin+box per sample, with stats in the tick labels and a Q10 line |

## Notes and pitfalls
- N50, median and coverage use **all** reads. Only the drawing is clipped at `--clip` (10 kb). The old script computed N50 on clipped reads, so its numbers are lower.
- Mean Q is the dorado `qs` tag. Without it, Q is computed from base qualities (error-probability mean), which can differ slightly from dorado because dorado skips the first 60 bases.
- Coverage uses the nuclear sacCer3 size (12.07 Mb), or `--genome-size`.
- Sample labels come from `--label-by auto`: the `BC` tag (e.g. `barcode09`), otherwise `barcodeNN` in the file name, otherwise the file stem.
- Secondary and supplementary records are skipped, so each read counts once.

## Provenance
Port of `process_nanopore_sequencing/Nov_20_25/no_trim/plot_read_data.py` (same styling). The input was the `samtools view | awk` `*_read_len_qual` files from `DM1746/script_sup_model.sh`. The filter cascade follows `bg_investigation/nov03/analyze_nov03.py`.
