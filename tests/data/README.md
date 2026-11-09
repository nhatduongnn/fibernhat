# Test data (small, committed on purpose)

| file | content |
|---|---|
| `nov20_barcode09_subset.bam` (3.8 MB) | Nov_20_25 run, barcode09 (dorado 0.9.1, sup@v5.0.0 + 6mA@v3, aligned to `ref/sacCer3/custom_sacCer3_lambda.fasta`). Contents: 30% of the chrI reads (seed 7), all reads at GAL10–GAL1 (chrII:277500-280000) and at TEA1 (chrXV:951269-957269), chrM:1-5000, and 1% of the unmapped reads |
| `nov20_barcode10_subset.bam` (5.6 MB) | the same for barcode10, with 15% of the chrI reads |
| `methylated_m6a_plotting_test.bam` (25 KB) | tiny fixture from `hmm_for_clustering_claude_code` |

These were made with `samtools view --subsample` from `process_nanopore_sequencing/Nov_20_25/no_trim/demux_output/*_barcode{09,10}.bam`, then `samtools cat`, `samtools sort`, and duplicates removed.
