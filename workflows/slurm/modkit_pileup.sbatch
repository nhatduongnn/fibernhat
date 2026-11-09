#!/bin/bash
# Standard modkit bedMethyl pileup (for sharing / other tools). fibernhat itself
# does not need it. Pinned modkit 0.4.4 (exact parity with fibernhat pysam engine).
#   sbatch workflows/slurm/modkit_pileup.sbatch BAM=x.bam OUT=x.bed [THR=0.8222656|none] [MODKIT=...]
#SBATCH --job-name=modkit_pileup
#SBATCH --partition=compsci
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=modkit_pileup_%j.out
set -euo pipefail
for kv in "$@"; do export "$kv"; done
MODKIT="${MODKIT:-${FIBERNHAT_MODKIT:-$(command -v modkit || echo /usr/xtmp/nd141/programs/dist_modkit_v0.4.4_7cf558c/modkit)}}"
: "${BAM:?}"; : "${OUT:?}"
if [ "${THR:-none}" = none ]; then FILT=(--no-filtering); else FILT=(--filter-threshold "A:$THR"); fi
"$MODKIT" pileup "$BAM" "$OUT" -t "${SLURM_CPUS_PER_TASK:-8}" "${FILT[@]}" --log-filepath "$OUT.log"
