#!/bin/bash
# Rebuild examples/gallery/*.png from the test BAMs (~2 min). Run from the repo root.
set -euo pipefail
F="${FIBERNHAT:-.venv/bin/fibernhat}"
B="tests/data/nov20_barcode09_subset.bam tests/data/nov20_barcode10_subset.bam"
T=$(mktemp -d); G=examples/gallery; mkdir -p $G
$F thresholds $B -o $T/thr.tsv
$F qc $B -o $T/qc --formats png
cp $T/qc/summary_read_length.png $G/qc_summary_read_length.png
cp $T/qc/summary_quality.png $G/qc_summary_quality.png
cp $T/qc/barcode10_hist.png $G/qc_hist_barcode10.png
$F fiberplot --bam $B --region chrII:277500-280000 --mark 278982 --title "GAL10-GAL1 (Nov_20_25 test subset)" \
   --thresholds $T/thr.tsv -o $T/fp --name GAL1 --formats png
cp $T/fp/GAL1.png $G/fiberplot_GAL1.png
$F fiberplot --bam $B --region chrXV:953300-955300 --mark 954213 --title "TEA1" --thresholds $T/thr.tsv \
   --sort meth -o $T/fp --name TEA1 --formats png
cp $T/fp/TEA1.png $G/fiberplot_TEA1.png
$F aggregate --bam $B --sites chereji --chrom chrI --half 500 --sort meth_in_win --thresholds $T/thr.tsv \
   -o $T/agg --name chereji_chrI --title "Chereji +1/-1 nucleosomes (chrI)" --formats png
cp $T/agg/chereji_chrI.png $G/aggregate_chereji_chrI.png
$F aggregate --bam $B --sites tss --chrom chrI --half 1000 --sort meth_in_win --thresholds $T/thr.tsv \
   -o $T/agg --name tss_chrI --title "Park TSS (chrI)" --formats png
cp $T/agg/tss_chrI.png $G/aggregate_tss_chrI.png
$F quintiles --bam $B --sites dyad,tss --thresholds $T/thr.tsv -o $T/q --formats png
cp $T/q/quintile_profiles.png $G/quintile_profiles.png
echo "gallery -> $G (tmp: $T)"
