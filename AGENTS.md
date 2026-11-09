# AGENTS.md: instructions for AI agents working in this repo

This repo is **fibernhat**, a set of tools for nanopore Fiber-seq (m6A) data. Everything reads BAM files directly; there is no parquet step and modkit is not required.

## Start here
1. Read `CAPABILITIES.md` to find the skill that matches the request.
2. Read that skill's `.claude/skills/<skill>/SKILL.md` and follow it. Each one gives the exact command, its parameters, and the common mistakes.
3. Run commands through the `fibernhat` CLI. Only import `fibernhat.*` in Python when the user needs a custom figure.

## Environment
- If `.venv/` exists, use `.venv/bin/fibernhat`.
- Otherwise, if `uv` is available: `uv venv .venv && uv pip install -e ".[all]"`.
- Otherwise use conda: `conda env create -f environment.yml && conda activate fibernhat`.
- Python 3.9 or newer. Dependencies: pysam, numpy, pandas, scipy, matplotlib, statsmodels, pyyaml.
- modkit is optional. Pin version 0.4.4, which is the version tested for exact parity. Set its path in `config/site.yaml` (`tools.modkit`) or in `FIBERNHAT_MODKIT`.

## Rules
- **Never commit raw data.** No BAM, pod5, parquet or large TSV goes in git. The only exception is the small BAMs in `tests/data/`. Write outputs to a directory outside the repo or under `outputs/`, which is gitignored.
- **Never commit secrets.** No GitHub tokens or credentials. Never copy another repo's `.git/`.
- **Use the shared core.** All MM/ML parsing, read filtering and m6A calling lives in `src/fibernhat/io/bam.py`. New code calls these functions (`iter_reads`, `read_to_fiber`, `ReadFilter`, `Thresholds`) instead of re-parsing tags.
- **m6A call default.** Each sample gets its own cutoff at the 10th percentile of call confidence (modkit `sample-probs` semantics). This estimate runs automatically, and `fibernhat thresholds` saves it to a TSV. Use `--threshold-ml N` only to reproduce an older analysis:
  - 192 for the paper fiberplots
  - 128 for modkit argmax and the original quintile analysis
- **Default read filters** (the lab "filtnordna" cascade): primary mapped reads, MAPQ ≥ 20, reference span ≥ 500, `qs` ≥ 10. chrM, chrLambda and rDNA are excluded. The CLI can change any of these.
- **Keep both engines in sync.** If you change anything that has both a pysam and a modkit engine (`thresholds`, `pileup`, `quintiles`), run `pytest tests/test_parity.py` (or `fibernhat parity --bam ...`). Results must be identical.
- **Coordinates.** Region strings are 1-based inclusive (`chrI:60500-63500`); code works in 0-based half-open coordinates. Chromosome names use UCSC style (`chrIV`). `refs.norm_chrom` converts `chr4`/`IV`/`4`.

## Adding a new skill
1. Copy `.claude/skills/_template/` to `.claude/skills/<name>/` and fill in `SKILL.md`.
2. Put the logic in `src/fibernhat/` and plotting in `src/fibernhat/plot/`, reusing `io/`.
3. Add a CLI subcommand in `src/fibernhat/cli.py`.
4. Add a row to `CAPABILITIES.md`; `tests/test_capabilities.py` enforces this.
5. Add a test in `tests/`, an example in `examples/`, and a line in `CHANGELOG.md`.
6. Record where the code came from (original script paths) in the SKILL.md under "Provenance".

## Layout
```
src/fibernhat/io/        bam.py (reader/filters/calls), thresholds.py, refs.py (site sets), mnase.py
src/fibernhat/qc/        read_stats.py
src/fibernhat/plot/      fibers.py, aggregate.py, metaprofile.py, qc.py, style.py
src/fibernhat/           aggregate.py, quintiles.py, pileup.py, cli.py, config.py
src/fibernhat/annotation SGD/Park gene track (vendored lab DMS_seq code)
ref/                     sacCer3+lambda FASTA, GFFs, site sets, lab threshold TSV (see ref/README.md)
workflows/slurm/         dorado basecall/align/demux, generic fibernhat job, modkit pileup
tests/                   pytest; tests/data has small Nov_20_25 subsets (chrI 15-30%, GAL1, TEA1 loci)
```
