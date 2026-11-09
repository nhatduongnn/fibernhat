---
name: basecall-align-dorado
description: Basecall nanopore pod5 with dorado sup + 6mA, align to sacCer3+lambda, demultiplex barcodes, then run fibernhat read QC and per-barcode thresholds on SLURM. Use when processing a new nanopore fiber-seq run from pod5.
---
# Dorado basecall → align → demux → QC

You don't need this if you already have BAMs: every other skill takes any aligned BAM that has MM/ML tags.

## Command
```bash
sbatch --gres=gpu:a5000:1 workflows/slurm/basecall_align_demux.sbatch \
   POD5_DIR=/path/to/pod5 OUTDIR=/path/to/out LABEL=DM1800 \
   KIT=SQK-RBK114-24 MODEL=sup,6mA DORADO=/path/to/dorado-2.0.1/bin/dorado
# optional: SOFTCLIP_SUPP=1 (minimap2 -Y), SKIP_BASECALL=1 (only QC an existing demux_output/)
```
The steps are:
1. `dorado basecaller sup,6mA --recursive --no-trim`
2. `dorado aligner <ref>`
3. `dorado demux --no-trim --sort-bam --kit-name <KIT>`
4. `samtools index`
5. `fibernhat qc` and `fibernhat thresholds --all-reads`

## Lab history (what the existing BAMs are)
| runs | dorado | models | notes |
|---|---|---|---|
| Mar_19_24 – Feb_25 | 0.5.3 | sup + 6mA, aligned during basecalling, kit SQK-NBD114-24 | `--min-qscore 9`, `demux --no-classify` |
| Mar_20_25 – Nov_20_25 | 0.9.1 | `sup@v5.0.0` + `6mA@v3` | `--no-trim`, separate aligner, RBK114-24, flat demux layout |
| Mar_09_26 | 1.4.0 | `sup@v5.2.0`, **no 6mA** | strain reference yAM095, `-Y` |
| DM1746 (Jun 2026) | 2.0.1 | `sup@v5.2.0` + `6mA@v1` | nested `bam_pass/barcodeNN/` demux layout |

## Notes
- Reference: `ref/sacCer3/custom_sacCer3_lambda.fasta`, which is what the lab BAMs were aligned to (`@SQ` md5 match). Dorado indexes it on the fly; there is no `.mmi`.
- Dorado ≥ 1.x writes demux output in nested folders. The script finds BAMs with `find`, so both layouts work.
- If the log says "Failed to find optimal batch size", add `--batchsize-benchmarks-file` to save GPU time.
- Pod5 files and BAMs must never be committed.

## Provenance
`process_nanopore_sequencing/DM1746/script_sup_model.sh` (current), `Nov_20_25/no_trim/script_sup_model.sh`, `Mar_09_26_*/no_trim/script_and_parquet_convert_final.sh` (`-Y`). The parquet conversion was intentionally dropped.
