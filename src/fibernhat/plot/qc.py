"""Read-length / quality plots. Port of process_nanopore_sequencing/Nov_20_25/
no_trim/plot_read_data.py (same layout and styling) with fixes:
N50/median on unclipped data, no stray 100/200 bp lines, genome size from
chrom.sizes, input straight from the BAM-derived read table."""
from __future__ import annotations

import os
from typing import List, Optional

import numpy as np
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D

from ..qc.read_stats import n50
from . import style

FONT_TITLE, FONT_LABEL, FONT_TICK, FONT_LEGEND = 18, 15, 13, 12


def _spines(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)


def _plt():
    return style.setup_mpl()


def _qc_style(fn):
    """Apply the original whitegrid look to QC figures only. plt.style.use()
    would change matplotlib globally and leak grids into later fiberplots."""
    import functools

    @functools.wraps(fn)
    def wrapper(*a, **k):
        import matplotlib.pyplot as plt
        with plt.style.context("seaborn-v0_8-whitegrid"):
            return fn(*a, **k)
    return wrapper


AUTO_BIN_MIN_READS = 20000


def auto_bin_width(n_reads: int) -> int:
    """10 bp bins (the original) for full runs; 100 bp when there are fewer
    than 20k reads. With sparse 10 bp bins a base-weighted histogram shows
    bars on the lines y = k * length (k = reads in the bin)."""
    return 10 if n_reads >= AUTO_BIN_MIN_READS else 100


@_qc_style
def plot_histogram(lengths, out_prefix, label, color, clip=10000, ymax=None,
                   log_x=False, formats=("png", "svg"), bin_width="auto"):
    plt = _plt()
    L = np.asarray(lengths, float)
    L = L[L > 0]
    N50, med = n50(L), float(np.median(L)) if L.size else 0
    arr = L[L < clip]
    fig, ax = plt.subplots(figsize=(9, 5))
    bw = auto_bin_width(L.size) if bin_width in (None, "auto") else int(bin_width)
    bins = np.arange(0, clip + bw, bw)
    ax.hist(arr, bins=bins, weights=arr, color=color, alpha=0.65, edgecolor="none", label="Basepairs")
    top = ymax or ax.get_ylim()[1]
    ax.axvline(N50, color="#333333", lw=1.6, ls="--")
    ax.axvline(med, color="#666666", lw=1.2, ls=":")
    ax.text(N50 + 100, 0.95 * top, f"N50 = {N50:,} bp", fontsize=FONT_LEGEND, va="top")
    ax.text(med + 100, 0.83 * top, f"Median = {med:,.0f} bp", fontsize=FONT_LEGEND, va="top")
    if ymax:
        ax.set_ylim(0, ymax)
    ax.set_xlim(0, clip)
    if log_x:
        ax.set_xscale("log"); ax.set_xlim(10, clip)
    else:
        ax.set_xticks(np.arange(0, clip + 1, 1000))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
        plt.setp(ax.get_xticklabels(), rotation=45)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{x / 1e6:.1f}M" if x >= 1e6 else f"{x / 1e3:.0f}k"))
    ax.set_xlabel("Read length (bp)", fontsize=FONT_LABEL)
    ax.set_ylabel("Total bases sequenced", fontsize=FONT_LABEL)
    ax.set_title(f"Bases sequenced — {label}\n(N50/median on all reads; {bw} bp bins, clipped at {clip:,} bp)",
                 fontsize=FONT_TITLE - 4)
    ax.tick_params(labelsize=FONT_TICK)
    ax.legend(loc="upper right", frameon=False, fontsize=FONT_LEGEND)
    _spines(ax)
    p = style.save(fig, out_prefix, formats, 300)
    plt.close(fig)
    return p


def _violin_box(ax, data, colors):
    pos = np.arange(1, len(data) + 1)
    vp = ax.violinplot(data, positions=pos, widths=0.7, showmeans=False,
                       showextrema=False, showmedians=False)
    for b, c in zip(vp["bodies"], colors):
        b.set_facecolor(c); b.set_edgecolor(c); b.set_alpha(0.35); b.set_linewidth(1.0)
    bp = ax.boxplot(data, positions=pos, widths=0.25, patch_artist=True, showfliers=False,
                    showmeans=True,
                    meanprops=dict(marker="D", markerfacecolor="white", markeredgecolor="black",
                                   markersize=6, zorder=5),
                    medianprops=dict(color="black", lw=2), whiskerprops=dict(lw=1.2),
                    capprops=dict(lw=1.5), boxprops=dict(lw=1.0))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.6)
    return pos


@_qc_style
def plot_summary(groups: List[tuple], out_dir: str, genome_size: int, clip=10000,
                 formats=("png", "svg")):
    """groups = [(label, lengths, quals), ...]"""
    plt = _plt()
    labels = [g[0] for g in groups]
    colors = [style.label_color(i) for i in range(len(groups))]
    n = len(groups)
    legend = [Line2D([0], [0], marker="D", color="w", markerfacecolor="white",
                     markeredgecolor="black", markersize=7, label="Mean"),
              Line2D([0], [0], color="black", lw=2, label="Median")]
    paths = []
    # length
    fig, ax = plt.subplots(figsize=(max(8, 2.5 * n), 7))
    data = [np.clip(np.asarray(g[1], float), 0, clip) for g in groups]
    pos = _violin_box(ax, data, colors)
    ticks = []
    for lab, L, _ in groups:
        L = np.asarray(L, float)
        ticks.append(f"{lab}\nn={len(L):,}\nN50={n50(L):,} bp\nMedian={np.median(L):,.0f} bp\n"
                     f"Mean={L.mean():,.0f} bp\nCov={L.sum() / genome_size:.1f}×")
    ax.set_xticks(pos); ax.set_xticklabels(ticks, fontsize=FONT_TICK - 2)
    ax.set_ylim(0, clip)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_ylabel("Read length (bp)", fontsize=FONT_LABEL)
    ax.set_title(f"Read length distribution (violin clipped at {clip:,} bp; stats unclipped)",
                 fontsize=FONT_TITLE - 2)
    ax.legend(handles=legend, loc="upper right", frameon=False, fontsize=FONT_LEGEND)
    _spines(ax)
    paths += style.save(fig, os.path.join(out_dir, "summary_read_length"), formats, 300)
    plt.close(fig)
    # quality
    fig, ax = plt.subplots(figsize=(max(8, 2.5 * n), 7))
    data = [np.asarray(g[2], float)[~np.isnan(np.asarray(g[2], float))] for g in groups]
    pos = _violin_box(ax, data, colors)
    ax.axhline(10, color="#e05252", lw=1.2, ls="--", alpha=0.7)
    ax.text(n + 0.55, 10.3, "Q10", color="#e05252", fontsize=FONT_LEGEND)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{lab}\nn={len(q):,}\nMedian Q={np.median(q):.1f}\nMean Q={np.mean(q):.1f}"
                        for lab, q in zip(labels, data)], fontsize=FONT_TICK - 2)
    ax.set_ylim(0, 50)
    ax.set_ylabel("Mean read quality (Phred)", fontsize=FONT_LABEL)
    ax.set_title("Read quality score distribution", fontsize=FONT_TITLE - 2)
    ax.legend(handles=legend, loc="upper right", frameon=False, fontsize=FONT_LEGEND)
    _spines(ax)
    paths += style.save(fig, os.path.join(out_dir, "summary_quality"), formats, 300)
    plt.close(fig)
    return paths


@_qc_style
def plot_cascade(casc, out_dir, formats=("png", "svg")):
    plt = _plt()
    labs = list(dict.fromkeys(casc.label))
    fig, ax = plt.subplots(figsize=(max(7, 1.6 * len(labs) + 3), 4.5))
    steps = [s for s in dict.fromkeys(casc.step) if s != "all_primary"]
    bottom = np.zeros(len(labs))
    tot = np.array([casc[(casc.label == l) & (casc.step == "all_primary")].remaining.iloc[0] for l in labs], float)
    cmap = plt.get_cmap("tab20")
    for i, s in enumerate(steps):
        v = np.array([casc[(casc.label == l) & (casc.step == s)].removed.sum() for l in labs]) / tot
        ax.bar(labs, v, bottom=bottom, color=cmap(i), label=s.replace("_", " "))
        bottom += v
    ax.bar(labs, 1 - bottom, bottom=bottom, color="#4c9a2a", label="PASS")
    ax.set_ylabel("Fraction of primary reads")
    ax.set_title("Read filter cascade")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    _spines(ax)
    p = style.save(fig, os.path.join(out_dir, "filter_cascade"), formats, 200)
    plt.close(fig)
    return p
