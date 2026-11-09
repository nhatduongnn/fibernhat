# Changelog

## 0.1.0 (2026-11-08)
- First release, BAM-only (pysam). No parquet.
- Skills: `read-qc`, `call-thresholds`, `raw-fiberplot`, `aggregate-fiberplot`, `quintile-metaprofile`, `basecall-align-dorado`.
- Default m6A call: a per-sample cutoff at the 10th percentile of call confidence (modkit `sample-probs` semantics).
- modkit is an optional engine, pinned to 0.4.4. `thresholds`, `pileup` and `quintiles` give identical results with the pysam and modkit engines on the test BAMs.
  - modkit 0.6.1 `pileup` leaves out a few calls that `modkit extract` reports (59 of 1.4M positions on the test data were 1 lower), so it is not used.
- Fixes compared with the old scripts:
  - Read QC: N50 is computed on unclipped reads.
  - Read QC: secondary and supplementary records are no longer counted.
  - Coverage uses the genome size from the chrom sizes file.
- Reference: sacCer3 plus lambda (matches the lab BAM `@SQ` lines), SGD R64-1-1 GFF, Park 2014 TSS/PAS, Brogaard/Chereji nucleosome sites, Rossi Abf1/Reb1 sites.
