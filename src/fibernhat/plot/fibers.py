"""Raw single-molecule fiberplots: gene track + (MNase) + m6A pileup + fibers.

Drawing conventions follow the paper figure modules
(hmm_for_clustering_claude_code/scripts/fig2_panels.py, fig1_panels.py):
one light-grey backbone per read, a steelblue dot per methylated A, strand
sidebar (+ navy / - dark red), MNase as a KDE-coloured V-plot, SGD ORFs with
Park 2014 TSS/PAS arrows (Dubious ORFs hidden).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..io.bam import Fiber
from . import style


# --------------------------------------------------------------------------
# ordering
# --------------------------------------------------------------------------
def order_fibers(fibers: List[Fiber], how: str = "start", seed: int = 0,
                 window: Optional[Tuple[int, int]] = None) -> List[int]:
    """Return row order (row 0 drawn at the TOP).

    how: start | shuffle | meth (most methylated first, within `window` if
    given) | length (longest first) | strand (+ then -, each by start)
    """
    n = len(fibers)
    if how == "shuffle":
        return list(np.random.default_rng(seed).permutation(n))
    if how == "meth":
        def frac(f):
            m = np.ones(f.pos.shape, bool) if window is None else (f.pos >= window[0]) & (f.pos < window[1])
            c = f.call[m]
            k = (c >= 0).sum()
            return (c == 1).sum() / k if k else -1
        return list(np.argsort([-frac(f) for f in fibers], kind="stable"))
    if how == "length":
        return list(np.argsort([-(f.end - f.start) for f in fibers], kind="stable"))
    if how == "strand":
        return sorted(range(n), key=lambda i: (fibers[i].strand != "+", fibers[i].start, fibers[i].end))
    return sorted(range(n), key=lambda i: (fibers[i].start, fibers[i].end))


# --------------------------------------------------------------------------
# panels
# --------------------------------------------------------------------------
def draw_gene_track(ax, chrom: str, start1: int, end1: int, fontsize: float = 7):
    """SGD ORFs + TSS/PAS for [start1, end1] (1-based). Returns False if the
    chromosome has no annotation."""
    from ..annotation import chrom_str_to_int, make_orf_annotation_plotter
    ci = chrom_str_to_int(chrom)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    if ci is None:
        ax.text(0.5, 0.5, f"no annotation for {chrom}", transform=ax.transAxes,
                ha="center", va="center", fontsize=6, color="grey")
        return False
    p = make_orf_annotation_plotter()
    p.set_span_chrom((start1, end1), int(ci))
    p.plot_orf_annotations(ax, orf_classes=["Verified", "Uncharacterized"])
    for t in ax.texts:
        t.set_fontsize(fontsize)
    ax.set_xlim(start1, end1)
    ax.set_xticks([]); ax.set_yticks([])
    return True


def draw_mnase(ax, mids, lens, start1, end1, max_len=300, lab_fs=style.LAB_FS,
               tick_fs=style.TICK_FS):
    style.plot_density_scatter(mids, lens, ax=ax, s=2, cmap="magma_r")
    ax.axhline(147, color="gray", lw=0.6, ls=":")
    ax.set_xlim(start1, end1)
    ax.set_ylim(0, max_len)
    ax.set_ylabel("MNase\nfrag (bp)", fontsize=lab_fs)
    ax.tick_params(labelsize=tick_fs, labelbottom=False)


def pileup_profile(fibers: List[Fiber], s0: int, e0: int, smooth: int = 31,
                   min_cov: int = 5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-bp m6A fraction over the region from the fibers' calls.
    Returns (x 1-based, frac (smoothed, NaN where coverage < min_cov), coverage)."""
    W = e0 - s0
    mod = np.zeros(W)
    val = np.zeros(W)
    for f in fibers:
        m = (f.pos >= s0) & (f.pos < e0) & (f.call >= 0)
        p = f.pos[m] - s0
        np.add.at(val, p, 1)
        np.add.at(mod, p[f.call[m] == 1], 1)
    k = np.ones(max(1, smooth))
    ms = np.convolve(mod, k, "same")
    vs = np.convolve(val, k, "same")
    with np.errstate(all="ignore"):
        frac = ms / vs
    # coverage in reads ~ called-A density; mask with a rough per-read count
    cov = np.zeros(W)
    for f in fibers:
        a, b = max(f.start, s0) - s0, min(f.end, e0) - s0
        if b > a:
            cov[a:b] += 1
    frac[cov < min_cov] = np.nan
    return np.arange(s0, e0) + 1, frac, cov


def draw_raw_fibers(ax, fibers: List[Fiber], order: Sequence[int], s0: int, e0: int,
                    dot_color=style.METH_COLOR, dot_s=None, line_lw=None):
    n = len(order)
    scale = min(1.0, 313.0 / max(n, 1))
    dot_s = dot_s if dot_s is not None else max(0.6, 4.0 * scale)
    line_lw = line_lw if line_lw is not None else max(0.08, 0.35 * scale)
    ys = np.arange(n)[::-1]  # row 0 at the top
    starts = np.array([fibers[i].start for i in order]) + 1
    ends = np.array([fibers[i].end for i in order])
    ax.hlines(ys, np.clip(starts, s0 + 1, e0), np.clip(ends, s0 + 1, e0),
              color=style.BACKBONE_COLOR, lw=line_lw, zorder=1)
    xs, yy = [], []
    for y, i in zip(ys, order):
        mp = fibers[i].meth_pos
        xs.append(mp + 1)
        yy.append(np.full(mp.size, y))
    if xs:
        ax.scatter(np.concatenate(xs), np.concatenate(yy), s=dot_s, c=dot_color,
                   marker=".", edgecolors="none", zorder=2, rasterized=True)
    ax.set_xlim(s0 + 1, e0)
    ax.set_ylim(-0.5, n - 0.5)
    ax.tick_params(axis="y", left=False, labelleft=False)


def draw_strand_sidebar(ax, fibers, order, ylabel, lab_fs=style.LAB_FS, tick_fs=style.TICK_FS):
    from matplotlib.colors import to_rgb
    n = len(order)
    rgb = np.array([[to_rgb(style.STRAND_COLOR.get(fibers[i].strand, "#888888"))] for i in order])
    ax.imshow(rgb, aspect="auto", interpolation="nearest", extent=[0, 1, n - 0.5, -0.5])
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([0, n - 1] if n > 1 else [0])
    ax.set_yticklabels(["1", str(n)] if n > 1 else ["1"])
    ax.set_ylabel(ylabel, fontsize=lab_fs, labelpad=1)
    ax.tick_params(labelsize=tick_fs, length=2, pad=1)
    for sp in ax.spines.values():
        sp.set_visible(False)


def draw_marks(ax, marks, color="black"):
    for x in marks or ():
        ax.axvline(x, color=color, ls="--", lw=0.6, alpha=0.8, zorder=4)


# --------------------------------------------------------------------------
# composite
# --------------------------------------------------------------------------
def fiberplot(fibers: List[Fiber], region: Tuple[str, int, int], out_prefix: str,
              title: str = "", sort: str = "start", seed: int = 0,
              mnase: Optional[pd.DataFrame] = None, marks: Sequence[int] = (),
              genes: bool = True, pileup: bool = True, split_by_label: bool = True,
              smooth: int = 31, width: float = 7.0, max_fiber_height: float = 7.0,
              formats=("png", "pdf"), dpi: int = 300, threshold_note: str = ""):
    """Draw and save the composite figure. Returns list of written paths."""
    plt = style.setup_mpl()
    from matplotlib.gridspec import GridSpec
    chrom, s0, e0 = region

    if split_by_label:
        labels = list(dict.fromkeys(f.label for f in fibers))
    else:
        labels = ["all"]
    groups = {lab: [f for f in fibers if (f.label == lab or not split_by_label)] for lab in labels}
    groups = {k: v for k, v in groups.items() if v}
    labels = list(groups)

    n_tot = max(1, sum(len(v) for v in groups.values()))
    fib_h = min(max_fiber_height, max(1.2, 0.012 * n_tot + 0.6 * len(labels)))
    rows, heights = [], []
    if genes:
        rows.append("genes"); heights.append(0.55)
    if mnase is not None:
        rows.append("mnase"); heights.append(0.9)
    if pileup:
        rows.append("pileup"); heights.append(0.8)
    for lab in labels:
        rows.append(("fib", lab))
        heights.append(max(0.35, fib_h * len(groups[lab]) / n_tot))
    fig = plt.figure(figsize=(width, sum(heights) + 0.7))
    gs = GridSpec(len(rows), 2, figure=fig, width_ratios=[0.025, 1], height_ratios=heights,
                  hspace=0.12, wspace=0.015)
    main_axes = []
    first_main = None
    for r, kind in enumerate(rows):
        ax = fig.add_subplot(gs[r, 1], sharex=first_main) if first_main else fig.add_subplot(gs[r, 1])
        if first_main is None and kind != "genes":
            first_main = ax
        if kind == "genes":
            draw_gene_track(ax, chrom, s0 + 1, e0)
            if title:
                ax.set_title(title, fontsize=style.TITLE_FS, loc="left")
        elif kind == "mnase":
            draw_mnase(ax, mnase["mid"].values, mnase["length"].values, s0 + 1, e0)
            draw_marks(ax, marks)
        elif kind == "pileup":
            for i, lab in enumerate(labels):
                x, fr, _ = pileup_profile(groups[lab], s0, e0, smooth=smooth)
                ax.plot(x, fr, lw=0.8, color=style.label_color(i), label=f"{lab} (n={len(groups[lab])})")
            ax.set_ylabel(f"m6A frac\n({smooth} bp)", fontsize=style.LAB_FS)
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=style.TICK_FS, labelbottom=False)
            ax.legend(fontsize=style.LEG_FS, frameon=False, loc="upper right", ncol=min(4, len(labels)))
            draw_marks(ax, marks)
            main_axes.append(ax)
        else:
            lab = kind[1]
            fl = groups[lab]
            order = order_fibers(fl, sort, seed, window=(s0, e0))
            draw_raw_fibers(ax, fl, order, s0, e0)
            draw_marks(ax, marks)
            side = fig.add_subplot(gs[r, 0])
            draw_strand_sidebar(side, fl, order, lab if split_by_label else "fibers")
            ax.tick_params(labelsize=style.TICK_FS, labelbottom=(r == len(rows) - 1))
            main_axes.append(ax)
    ax = main_axes[-1] if main_axes else first_main
    ax.set_xlabel(f"{chrom} position (bp)", fontsize=style.LAB_FS)
    ax.set_xlim(s0 + 1, e0)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    note = f"m6A call: {threshold_note}" if threshold_note else ""
    if note:
        fig.text(0.99, 0.003, note, ha="right", va="bottom", fontsize=5, color="grey")
    paths = style.save(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return paths


def fiber_table(fibers: List[Fiber]) -> pd.DataFrame:
    return pd.DataFrame([dict(read_id=f.read_id, label=f.label, chrom=f.chrom, start=f.start,
                              end=f.end, strand=f.strand, mapq=f.mapq, n_meth=f.n_meth,
                              n_called=f.n_called, rate=f.rate) for f in fibers])
