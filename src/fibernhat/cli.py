"""fibernhat command line: `fibernhat <command> -h` for details.

Commands
  qc          read length / quality / filter-cascade plots + tables from BAM
  thresholds  per-sample 10th-percentile-confidence m6A call thresholds
  fiberplot   raw single-molecule fiberplot of a region (+ genes, MNase, pileup)
  aggregate   piled-up fibers + m6A profile around a site set (dyad/tss/abf1/...)
  quintiles   per-read m6A-rate quintiles -> metaprofiles at dyad/TSS/TF sites
  pileup      per-position m6A counts (bedMethyl-like) -- pysam or modkit
  parity      check the pysam engine against modkit on your BAMs
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)


# --------------------------------------------------------------------------
# shared options
# --------------------------------------------------------------------------
def _none_or_float(v):
    return None if str(v).lower() in ("none", "off", "") else float(v)


def add_filter_args(p):
    g = p.add_argument_group("read filters (defaults = lab 'filtnordna' cascade)")
    g.add_argument("--min-mapq", type=int, default=None, help="default 20")
    g.add_argument("--min-span", type=int, default=None, help="min reference span, default 500")
    g.add_argument("--min-qs", type=str, default=None,
                   help="min dorado mean Q (qs tag), default 10; 'none' disables")
    g.add_argument("--keep-chroms", action="store_true", help="do not exclude chrM / chrLambda")
    g.add_argument("--keep-rdna", action="store_true", help="do not exclude rDNA chrXII:451418-468929")
    g.add_argument("--keep-reads", help="file with read IDs to keep (one per line)")
    g.add_argument("--no-filter", action="store_true", help="primary mapped reads only, no other filter")


def add_label_args(p):
    p.add_argument("--label-by", default=None, choices=["auto", "file", "bc", "rg", "rg+bc", "none"],
                   help="how reads are grouped into samples (default auto: BC tag barcode, else "
                        "barcodeNN from file name, else file stem)")


def add_threshold_args(p):
    g = p.add_argument_group("m6A call threshold (default: per-sample 10th-pct confidence, estimated on the fly)")
    g.add_argument("--thresholds", help="TSV from `fibernhat thresholds` (or the lab roboNhat TSV)")
    g.add_argument("--threshold-ml", type=int, help="fixed binary call ML >= N (e.g. 192 = paper fiberplots, "
                                                    "128 = modkit argmax / original quintile analysis)")


def build_filter(a, cfg):
    from .io.bam import ReadFilter
    if getattr(a, "no_filter", False):
        return ReadFilter.none()
    f = ReadFilter.from_dict(cfg.get("filters"))
    if a.min_mapq is not None:
        f.min_mapq = a.min_mapq
    if a.min_span is not None:
        f.min_span = a.min_span
    if a.min_qs is not None:
        f.min_qs = _none_or_float(a.min_qs)
    if a.keep_chroms:
        f.exclude_chroms = ()
    if a.keep_rdna:
        f.exclude_regions = ()
    if a.keep_reads:
        with open(a.keep_reads) as fh:
            f.keep_read_ids = {l.strip() for l in fh if l.strip()}
    return f


def resolve_thresholds(a, bams, filt, label_by, cfg):
    from .io import thresholds as T
    th = T.resolve(a.thresholds, a.threshold_ml, bams, filt, label_by,
                   max_reads_per_label=cfg.get("threshold_estimate_reads", 20000))
    note = th.source if a.threshold_ml is not None else (
        th.source + ": " + ", ".join(f"{k} ML>={v.ml_meth}" for k, v in list(th.by_label.items())[:6]))
    print(f"[fibernhat] m6A call: {note}", file=sys.stderr)
    return th, note


def _out(a, default):
    os.makedirs(a.outdir or default, exist_ok=True)
    return a.outdir or default


def _mnase_path(a, cfg):
    if getattr(a, "no_mnase", False):
        return None
    p = a.mnase_bam or cfg.get("mnase_bam")
    if p and not os.path.exists(p):
        print(f"[fibernhat] MNase BAM not found: {p} (skipping track)", file=sys.stderr)
        return None
    return p


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------
def cmd_qc(a, cfg):
    from .io import refs
    from .qc import read_stats as R
    from .plot import qc as PQ
    out = _out(a, "fibernhat_qc")
    filt = build_filter(a, cfg)
    df = R.read_table(a.bam, a.label_by or cfg["label_by"], filt, a.max_reads)
    print(f"[fibernhat] {len(a.bam)} BAM(s), samples: {', '.join(df.label.unique())}", file=sys.stderr)
    gsize = a.genome_size or refs.nuclear_genome_size()
    df.to_csv(os.path.join(out, "read_stats.tsv.gz"), sep="\t", index=False, float_format="%.3f")
    summ = R.summary(df, gsize)
    summ.to_csv(os.path.join(out, "summary.tsv"), sep="\t", index=False, float_format="%.4f")
    casc = R.cascade(df)
    casc.to_csv(os.path.join(out, "filter_cascade.tsv"), sep="\t", index=False)
    fm = tuple(a.formats.split(","))
    groups = []
    for i, (lab, g) in enumerate(df.groupby("label", sort=False)):
        g2 = g if not a.pass_only else g[g.fail_reason == ""]
        PQ.plot_histogram(g2.read_len.values, os.path.join(out, f"{lab}_hist"), lab,
                          PQ.style.label_color(i), clip=a.clip, ymax=a.hist_ymax, log_x=a.log_x, formats=fm)
        groups.append((lab, g2.read_len.values, g2.mean_q.values))
    PQ.plot_summary(groups, out, gsize, clip=a.clip, formats=fm)
    PQ.plot_cascade(casc, out, formats=fm)
    print(summ.to_string(index=False))
    print(f"\nwrote {out}/")


def cmd_thresholds(a, cfg):
    from .io import thresholds as T
    filt = build_filter(a, cfg)
    df = T.estimate(a.bam, filt, a.label_by or cfg["label_by"], a.percentile,
                    None if a.all_reads else a.max_reads, engine=a.engine,
                    modkit=a.modkit or cfg.get("tools", {}).get("modkit"), threads=a.threads)
    T.write_tsv(df, a.out, note=f"engine={a.engine} bams={','.join(map(os.path.basename, a.bam))}")
    print(df.to_string(index=False))
    print(f"\nwrote {a.out}")


def cmd_fiberplot(a, cfg):
    from .io.bam import fetch_fibers, parse_region
    from .io import mnase as M
    from .plot.fibers import fiberplot, fiber_table
    out = _out(a, "fibernhat_fiberplot")
    filt = build_filter(a, cfg)
    label_by = a.label_by or cfg["label_by"]
    th, note = resolve_thresholds(a, a.bam, filt, label_by, cfg)
    marks = [int(x.replace(",", "")) for x in a.mark.split(",")] if a.mark else []
    if a.gene:
        from .annotation import find_gene
        g = find_gene(a.gene)
        a.region = f"{g['chrom']}:{max(1, g['tss'] - a.flank)}-{g['tss'] + a.flank}"
        if not a.mark:
            marks = [g["tss"]]
        a.name = a.name or g["name"]
        print(f"[fibernhat] {g['name']}: {a.region}, TSS {g['tss']} ({g['tss_source']})", file=sys.stderr)
    if not a.region:
        sys.exit("give --region or --gene")
    region = parse_region(a.region)
    fibers = fetch_fibers(a.bam, region, th, filt, label_by, require_span=a.require_span,
                          max_reads=a.max_reads, seed=a.seed)
    if not fibers:
        sys.exit(f"no reads pass the filters in {a.region}")
    mp = _mnase_path(a, cfg)
    mn = M.fragments(mp, *region, max_len=300) if mp else None
    name = a.name or a.region.replace(":", "_").replace(",", "")
    panels = [x for x in a.panels.split(",")
              if not ((x == "genes" and a.no_genes) or (x == "mnase" and a.no_mnase)
                      or (x == "pileup" and a.no_pileup))]
    paths = fiberplot(fibers, region, os.path.join(out, name), title=a.title or "",
                      sort=a.sort, seed=a.seed, mnase=mn, marks=marks, panels=panels,
                      split_by_label=not a.pool, smooth=a.smooth, width=a.width,
                      formats=tuple(a.formats.split(",")), threshold_note=note,
                      show_note=a.call_note)
    fiber_table(fibers).to_csv(os.path.join(out, name + "_fibers.tsv"), sep="\t", index=False,
                               float_format="%.4f")
    print("\n".join(paths))


def cmd_aggregate(a, cfg):
    from .io import refs, mnase as M
    from . import aggregate as A
    from .plot.aggregate import aggregate_plot
    out = _out(a, "fibernhat_aggregate")
    filt = build_filter(a, cfg)
    label_by = a.label_by or cfg["label_by"]
    th, note = resolve_thresholds(a, a.bam, filt, label_by, cfg)
    sites = refs.load_sites(a.sites)
    if a.chrom:
        sites = sites[sites.chrom.isin(a.chrom.split(","))]
    if a.max_sites and len(sites) > a.max_sites:
        sites = sites.sample(a.max_sites, random_state=a.seed)
    if a.no_orient:
        sites = sites.assign(sign=1)
    half = a.half or refs.DEFAULT_HALF.get(a.sites, 500)
    sf = A.collect(a.bam, sites, half, th, filt, label_by, require_span=not a.no_require_span,
                   max_per_site=a.max_per_site, seed=a.seed)
    if not len(sf.rows):
        sys.exit("no fibers overlap the sites")
    name = a.name or refs.site_label(a.sites)
    sf.save(os.path.join(out, name + "_fibers.npz"))
    A.profile_table(sf, name).to_csv(os.path.join(out, name + "_profile.tsv"), sep="\t", index=False,
                                     float_format="%.5f")
    mp = _mnase_path(a, cfg)
    mn = M.fragments_around_sites(mp, sites, half, max_len=250) if mp else None
    paths = aggregate_plot(sf, os.path.join(out, name), title=a.title or "", sort=a.sort,
                           max_rows=a.max_rows, seed=a.seed, smooth_bp=a.smooth, mnase=mn,
                           split_by_label=not a.pool, panels=a.panels,
                           formats=tuple(a.formats.split(",")),
                           xlabel=f"distance to {name} (bp)", threshold_note=note,
                           show_note=a.call_note)
    print("\n".join(paths))


def cmd_quintiles(a, cfg):
    from . import quintiles as Q
    from .plot.metaprofile import plot_quintile_profiles
    out = _out(a, "fibernhat_quintiles")
    filt = build_filter(a, cfg)
    label_by = a.label_by or cfg["label_by"]
    th, note = resolve_thresholds(a, a.bam, filt, label_by, cfg)
    sites = a.sites.split(",")
    halves = {}
    for kv in (a.half or "").split(","):
        if "=" in kv:
            k, v = kv.split("=")
            halves[k] = int(v)
    reads, prof, summ = Q.run(a.bam, th, filt, sites, halves, label_by, a.min_len, a.engine,
                              a.modkit or cfg.get("tools", {}).get("modkit"), a.threads)
    reads.to_csv(os.path.join(out, "per_read_rates.tsv.gz"), sep="\t", index=False, float_format="%.5f")
    prof.to_csv(os.path.join(out, "metaprofiles_raw.tsv.gz"), sep="\t", index=False, float_format="%.6g")
    summ.to_csv(os.path.join(out, "quintile_summary.tsv"), sep="\t", index=False, float_format="%.6g")
    Q.scores(prof, summ).to_csv(os.path.join(out, "footprint_scores.tsv"), sep="\t", index=False,
                                float_format="%.5g")
    paths = plot_quintile_profiles(prof, summ, os.path.join(out, "quintile_profiles"),
                                   smooth_bp=a.smooth, formats=tuple(a.formats.split(",")),
                                   threshold_note=note, show_note=a.call_note)
    print(summ.to_string(index=False))
    print("\n".join(paths))


def cmd_pileup(a, cfg):
    from . import pileup as P
    from .io.bam import parse_region
    filt = build_filter(a, cfg)
    label_by = a.label_by or cfg["label_by"]
    th, _ = resolve_thresholds(a, a.bam, filt, label_by, cfg)
    arr = P.pileup(a.bam, th, filt, label_by, engine=a.engine,
                   modkit=a.modkit or cfg.get("tools", {}).get("modkit"), threads=a.threads,
                   region=parse_region(a.region) if a.region else None)
    P.write_bedmethyl_like(arr, a.out)
    print(f"wrote {a.out}  (chrom, start, end, strand, n_valid, n_mod, percent)")


def cmd_parity(a, cfg):
    """pysam vs modkit: thresholds and per-position pileup must be identical."""
    from . import pileup as P
    from .io import thresholds as T
    from .io.bam import Thresholds
    mk = a.modkit or cfg.get("tools", {}).get("modkit")
    filt = build_filter(a, cfg)
    label_by = a.label_by or cfg["label_by"]
    rep = {}
    t1 = T.estimate(a.bam, filt, label_by, max_reads_per_label=None)
    t2 = T.estimate(a.bam, filt, label_by, max_reads_per_label=None, engine="modkit", modkit=mk)
    m = t1.merge(t2, on="label", suffixes=("_pysam", "_modkit"))
    rep["thresholds_identical"] = bool((m.ml_meth_pysam == m.ml_meth_modkit).all())
    print(m[["label", "ml_meth_pysam", "ml_meth_modkit", "n_reads_pysam"]].to_string(index=False))
    from .io.bam import CallThreshold
    th = Thresholds({r.label: CallThreshold(int(r.ml_meth), int(r.ml_canon)) for r in t1.itertuples()})
    for name, t in (("per_label", th), ("fixed_ml128", Thresholds.fixed(128))):
        d = P.compare(P.pileup(a.bam, t, filt, label_by),
                      P.pileup(a.bam, t, filt, label_by, engine="modkit", modkit=mk))
        rep[f"pileup_{name}_mismatches"] = int(len(d))
        if len(d):
            d.to_csv(f"parity_mismatch_{name}.tsv", sep="\t", index=False)
    print(json.dumps(rep, indent=2))
    ok = rep["thresholds_identical"] and all(v == 0 for k, v in rep.items() if k.endswith("mismatches"))
    print("PARITY OK" if ok else "PARITY FAILED")
    sys.exit(0 if ok else 1)


# --------------------------------------------------------------------------
def main(argv=None):
    from .config import load_config
    ap = argparse.ArgumentParser(prog="fibernhat", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", help="extra YAML config (merged over config/defaults.yaml + config/site.yaml)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("qc", help="read length/quality QC from BAM")
    p.add_argument("bam", nargs="+")
    p.add_argument("-o", "--outdir")
    p.add_argument("--clip", type=int, default=10000, help="x-limit for histograms / violins (stats unclipped)")
    p.add_argument("--hist-ymax", type=float)
    p.add_argument("--log-x", action="store_true")
    p.add_argument("--genome-size", type=int, help="for coverage; default nuclear sacCer3 size")
    p.add_argument("--pass-only", action="store_true", help="plot only reads passing the filters")
    p.add_argument("--max-reads", type=int, help="per BAM, for a quick look")
    p.add_argument("--formats", default="png,svg")
    add_filter_args(p); add_label_args(p)
    p.set_defaults(func=cmd_qc)

    p = sub.add_parser("thresholds", help="per-sample m6A call thresholds")
    p.add_argument("bam", nargs="+")
    p.add_argument("-o", "--out", default="m6a_thresholds.tsv")
    p.add_argument("--percentile", type=float, default=0.10)
    p.add_argument("--max-reads", type=int, default=50000, help="reads per sample (default 50000)")
    p.add_argument("--all-reads", action="store_true", help="use every read (as the lab TSV did)")
    p.add_argument("--engine", choices=["pysam", "modkit"], default="pysam")
    p.add_argument("--modkit"); p.add_argument("--threads", type=int, default=4)
    add_filter_args(p); add_label_args(p)
    p.set_defaults(func=cmd_thresholds)

    p = sub.add_parser("fiberplot", help="raw single-molecule fiberplot of a region")
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("--region", help="chrI:60500-63500 (1-based, inclusive)")
    p.add_argument("--gene", help="gene name or ORF (e.g. GAL1, YBR020W) instead of --region; the TSS is "
                                  "marked automatically (Park 2014 TSS, else the start codon)")
    p.add_argument("--flank", type=int, default=1000, help="bp around the TSS for --gene (default 1000)")
    p.add_argument("-o", "--outdir"); p.add_argument("--name"); p.add_argument("--title")
    p.add_argument("--sort", default="shuffle", choices=["shuffle", "start", "meth", "length", "strand"],
                   help="fiber row order (default shuffle, seed --seed, as in the lab figures)")
    p.add_argument("--panels", default="genes,mnase,fibers,pileup",
                   help="panels top to bottom, any subset/order of genes,mnase,fibers,pileup "
                        "(default genes,mnase,fibers,pileup; mnase skipped without an MNase BAM)")
    p.add_argument("--max-reads", type=int, help="random subsample (per region, all samples)")
    p.add_argument("--require-span", action="store_true", help="only reads spanning the whole region")
    p.add_argument("--pool", action="store_true", help="one fiber panel for all samples")
    p.add_argument("--mark", help="comma-separated positions to mark (dashed lines), e.g. a TSS")
    p.add_argument("--mnase-bam", help="paired-end MNase BAM (default: config mnase_bam)")
    p.add_argument("--no-mnase", action="store_true"); p.add_argument("--no-genes", action="store_true")
    p.add_argument("--no-pileup", action="store_true")
    p.add_argument("--smooth", type=int, default=31, help="pileup box smoothing (bp)")
    p.add_argument("--width", type=float, default=7.0)
    p.add_argument("--seed", type=int, default=0); p.add_argument("--formats", default="png,pdf")
    p.add_argument("--call-note", action="store_true", help="print the m6A cutoff used in a small footer")
    add_filter_args(p); add_label_args(p); add_threshold_args(p)
    p.set_defaults(func=cmd_fiberplot)

    p = sub.add_parser("aggregate", help="aggregate fiberplot around a site set")
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("--sites", required=True,
                   help="dyad | tss | abf1 | reb1 | chereji | tf | tf:<TF> | path (TSV/BED)")
    p.add_argument("--half", type=int, help="half window (bp); default per site set")
    p.add_argument("-o", "--outdir"); p.add_argument("--name"); p.add_argument("--title")
    p.add_argument("--sort", default="meth_in_win", choices=["meth_in_win", "shuffle", "meth_frac", "site", "none"],
                   help="fiber row order (default meth_in_win: most methylated within +/-100 bp first)")
    p.add_argument("--panels", default="mnase,fibers,counts,ratio",
                   help="panels top to bottom, any subset/order of mnase,fibers,counts,ratio "
                        "(counts = methylation pileup, ratio = methylated As / possible As)")
    p.add_argument("--max-rows", type=int, default=3000, help="fibers drawn per sample (profile uses all)")
    p.add_argument("--max-per-site", type=int, help="cap fibers per site")
    p.add_argument("--max-sites", type=int, help="random subset of sites")
    p.add_argument("--chrom", help="restrict to chromosomes (comma-separated)")
    p.add_argument("--no-require-span", action="store_true",
                   help="include reads that only partly cover a window")
    p.add_argument("--no-orient", action="store_true", help="do not flip minus-strand sites")
    p.add_argument("--pool", action="store_true")
    p.add_argument("--mnase-bam"); p.add_argument("--no-mnase", action="store_true")
    p.add_argument("--smooth", type=int, default=1, help="rolling mean (bp) on counts/ratio; 1 = raw (original)")
    p.add_argument("--seed", type=int, default=0); p.add_argument("--formats", default="png,pdf")
    p.add_argument("--call-note", action="store_true", help="print the m6A cutoff used in a small footer")
    add_filter_args(p); add_label_args(p); add_threshold_args(p)
    p.set_defaults(func=cmd_aggregate)

    p = sub.add_parser("quintiles", help="per-read m6A-rate quintile metaprofiles")
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("--sites", default="dyad,tss,abf1,reb1")
    p.add_argument("--half", help="override half windows, e.g. dyad=1000,tss=1500")
    p.add_argument("--min-len", type=int, default=500, help="min read length for quintile assignment")
    p.add_argument("-o", "--outdir")
    p.add_argument("--engine", choices=["pysam", "modkit"], default="pysam")
    p.add_argument("--modkit"); p.add_argument("--threads", type=int, default=4)
    p.add_argument("--smooth", type=int, default=21)
    p.add_argument("--formats", default="png,pdf")
    p.add_argument("--call-note", action="store_true", help="print the m6A cutoff used in a small footer")
    add_filter_args(p); add_label_args(p); add_threshold_args(p)
    p.set_defaults(func=cmd_quintiles)

    p = sub.add_parser("pileup", help="per-position m6A counts")
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--region")
    p.add_argument("--engine", choices=["pysam", "modkit"], default="pysam")
    p.add_argument("--modkit"); p.add_argument("--threads", type=int, default=4)
    add_filter_args(p); add_label_args(p); add_threshold_args(p)
    p.set_defaults(func=cmd_pileup)

    p = sub.add_parser("parity", help="verify pysam == modkit on these BAMs")
    p.add_argument("--bam", nargs="+", required=True)
    p.add_argument("--modkit")
    add_filter_args(p); add_label_args(p)
    p.set_defaults(func=cmd_parity)

    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    a.func(a, cfg)


if __name__ == "__main__":
    main()
