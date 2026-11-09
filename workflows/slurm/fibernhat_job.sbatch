#!/bin/bash
# Run any fibernhat command on a CPU node.
#   sbatch workflows/slurm/fibernhat_job.sbatch quintiles --bam a.bam b.bam -o out/
#SBATCH --job-name=fibernhat
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=fibernhat_%j.out
set -euo pipefail
REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
FIBERNHAT="${FIBERNHAT:-$REPO/.venv/bin/fibernhat}"
command -v "$FIBERNHAT" >/dev/null 2>&1 || FIBERNHAT=fibernhat
export MPLBACKEND=Agg
echo "[$(date)] fibernhat $*"
"$FIBERNHAT" "$@"
echo "[$(date)] done"
