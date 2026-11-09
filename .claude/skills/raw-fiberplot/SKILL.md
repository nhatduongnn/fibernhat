---
name: raw-fiberplot
description: Plot raw single-molecule fiber-seq (m6A) fiberplots of a genomic region straight from BAM, with SGD gene track, optional MNase V-plot, per-sample m6A pileup and per-read dot rows. Use for "fiberplot", "single-molecule plot", "show reads at GAL1 / a locus", locus views.
---
# Raw fiberplot (BAM → locus figure)

## Command
```bash
# by gene: window = TSS ± --flank (1000), TSS marked automatically
.venv/bin/fibernhat fiberplot --bam a.bam b.bam --gene TEA1 -o outputs/fiberplots
# by region:
.venv/bin/fibernhat fiberplot --bam a.bam b.bam --region chrII:277500-280000 \
    --mark 278982 --title "GAL10-GAL1" -o outputs/fiberplots
# options: --panels genes,mnase,fibers,pileup (any subset, any order, top to bottom)
#          --sort shuffle|start|meth|length|strand  --max-reads 300  --require-span
#          --pool (one fiber panel)  --no-mnase  --no-genes  --no-pileup  --smooth 31  --call-note
#          --threshold-ml 192 (paper figures)  --thresholds thr.tsv
```

## Figure, top to bottom (default `--panels genes,mnase,fibers,pileup`, the lab layout)
1. **Gene track:** SGD R64-1-1 ORFs (Verified and Uncharacterized; Dubious hidden) with Park 2014 TSS/PAS arrows. Watson genes are blue, Crick genes red. No box, ylim (−60, 118), as in `fig2_panels`/`build_tea1_plots.py`.
2. **MNase V-plot** (only if `mnase_bam` is set in `config/site.yaml` or `--mnase-bam` is given): fragment midpoint vs length, coloured by 2-D KDE, with a dotted line at 147 bp. Y axis on the right (0/150/300).
3. **Fibers:** one panel per sample (`--pool` for one panel, as in the paper figures). Rows are shuffled with a fixed seed. Each read is a grey line with a steelblue dot per methylated A. The left sidebar shows strand (+ navy, − dark red).
4. **m6A pileup:** per sample, methylated/called A per bp with a 31 bp box smoothing. Positions covered by fewer than 5 reads are hidden. This panel is a fibernhat addition (not in the paper locus figures); drop it with `--no-pileup`.

No title unless `--title` is given. The m6A cutoff used is printed in the terminal; `--call-note` also writes it in a small footer.

The command also writes `<name>_fibers.tsv` with per-read start/end/strand/n_meth/n_called/rate.

## Notes
- **TSS used by `--gene`:** the Park 2014 TSS when there is one, otherwise the start codon. The title says which. Some genes, including **GAL1 and GAL10, have no Park TSS** (they are silent in glucose). The lab's GAL1 TSS is **chrII:278982**, 39 bp upstream of the ATG at 279021, from `hmm_for_clustering_claude_code/docs/gal1_induction_clustering_handoff.md`. Use `--region chrII:277500-280000 --mark 278982` for the paper GAL10–GAL1 view.
- Useful loci from the lab figures:
  - GAL10–GAL1: `chrII:277500-280000`, TSS 278982, UAS window 278500–278750
  - ERV46: `chrI:60500-63500`, Abf1 sites 61170 and 62664
  - TEA1: `chrXV:951269-957269`, TSS 954213
- The m6A call defaults to a per-sample 10th-percentile-confidence cutoff, estimated on the fly from 20k reads per sample. For consistent figures across runs, make the TSV once with `fibernhat thresholds` and pass `--thresholds`.
- Reverse-strand reads are scored at reference T (template A), the same convention as modkit.
- Only methylated calls are drawn. No-calls and canonical A's are not drawn.
- `--region` is 1-based inclusive. Very deep regions get slow to draw; use `--max-reads`.
- The paper figures used ML ≥ 192, a seeded shuffle and one pooled panel without the pileup: `--threshold-ml 192 --pool --no-pileup`.

## Python
```python
from fibernhat.io.bam import fetch_fibers, parse_region, ReadFilter
from fibernhat.io.thresholds import resolve
from fibernhat.plot.fibers import fiberplot
```

## Provenance
The drawing code is a port of `hmm_for_clustering_claude_code/scripts/fig2_panels.py` (`draw_raw_fibers`, `draw_gene_track`, `draw_mnase_scatter`, `draw_sidebar`, `draw_locus_overlays`) and `fig1_panels.py` (`draw_panelG`). The same code is used in `collab/Strahl_lab/scripts/build_tea1_plots.py` and roboNhat `singlefiber/fig_permol_fibers.py`. The gene track is vendored from `DMS_seq/src/script/{orf_annotation_plotter,reference_data}.py`. The KDE scatter comes from DMS_seq `plot_utils.plot_density_scatter`, and the MNase reader from DMS_seq `read_bam.read_mnase_bam_single_chrom`. The BAM reader replaces `hmm .../src/bam_processor.py`, and its calls follow roboNhat `singlefiber/read_calls.py`.
