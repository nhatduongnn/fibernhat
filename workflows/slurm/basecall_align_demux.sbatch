#!/bin/bash
# Dorado basecall (sup + 6mA) -> align -> demux -> fibernhat read QC.
# Parameterized from process_nanopore_sequencing/DM1746/script_sup_model.sh (dorado 2.0.1)
# + the --mm2-opts "-Y" option of Mar_09_26/no_trim/script_and_parquet_convert_final.sh.
#
# Usage (edit the SBATCH header for your cluster, or override on the command line):
#   sbatch --gres=gpu:a5000:1 workflows/slurm/basecall_align_demux.sbatch \
#       POD5_DIR=/path/to/pod5 OUTDIR=/path/to/run_out LABEL=myrun \
#       [KIT=SQK-RBK114-24] [MODEL=sup,6mA] [REF=ref/sacCer3/custom_sacCer3_lambda.fasta] \
#       [DORADO=/path/to/dorado] [SOFTCLIP_SUPP=1] [SKIP_BASECALL=1] [SKIP_QC=1]
#
#SBATCH --job-name=fibernhat_dorado
#SBATCH --gres=gpu:a5000:1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --output=fibernhat_dorado_%j.out
#SBATCH --error=fibernhat_dorado_%j.err
set -euo pipefail
for kv in "$@"; do export "$kv"; done

REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)}"
[ -f "$REPO/pyproject.toml" ] || REPO="${SLURM_SUBMIT_DIR:-$PWD}"
: "${POD5_DIR:?set POD5_DIR=...}"
: "${OUTDIR:?set OUTDIR=...}"
LABEL="${LABEL:-run}"
KIT="${KIT:-SQK-RBK114-24}"
MODEL="${MODEL:-sup,6mA}"
REF="${REF:-$REPO/ref/sacCer3/custom_sacCer3_lambda.fasta}"
DORADO="${DORADO:-$(command -v dorado || echo /usr/xtmp/nd141/programs/dorado/dorado-2.0.1-linux-x64/bin/dorado)}"
FIBERNHAT="${FIBERNHAT:-$REPO/.venv/bin/fibernhat}"
mkdir -p "$OUTDIR"; cd "$OUTDIR"

echo "[$(date)] dorado: $($DORADO --version 2>&1 | head -1)  model: $MODEL  kit: $KIT"
if [ -z "${SKIP_BASECALL:-}" ]; then
  $DORADO basecaller "$MODEL" "$POD5_DIR" --recursive --no-trim > "${LABEL}_no_align.bam"
  MM2=(); [ -n "${SOFTCLIP_SUPP:-}" ] && MM2=(--mm2-opts "-Y")
  $DORADO aligner "${MM2[@]}" "$REF" "${LABEL}_no_align.bam" > "${LABEL}_aligned.bam"
  $DORADO demux --no-trim --sort-bam --kit-name "$KIT" --output-dir demux_output "${LABEL}_aligned.bam"
fi

# Per-barcode BAMs: dorado <=0.9 writes flat demux_output/*_barcodeNN.bam,
# dorado >=1.x writes demux_output/<exp>/<sample>/<run>/bam_pass/barcodeNN/*.bam
mapfile -t BAMS < <(find demux_output -name "*.bam" ! -name "*unclassified*" | sort)
echo "found ${#BAMS[@]} demuxed BAM(s)"
for b in "${BAMS[@]}"; do [ -f "$b.bai" ] || samtools index "$b"; done

if [ -z "${SKIP_QC:-}" ] && [ ${#BAMS[@]} -gt 0 ]; then
  "$FIBERNHAT" qc "${BAMS[@]}" -o "qc_${LABEL}"
  "$FIBERNHAT" thresholds "${BAMS[@]}" -o "${LABEL}_m6a_thresholds.tsv" --all-reads
fi
echo "[$(date)] done -> $OUTDIR"
