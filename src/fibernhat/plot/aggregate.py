"""Aggregate fiberplot: (MNase V-plot) + piled-up fibers + m6A fraction profile."""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from ..aggregate import SiteFibers, order_rows
from . import style


def smooth(y, k):
    if k <= 1:
        return y
    return pd.Series(y).rolling(k, center=True, min_periods=1).mean().values


def aggregate_plot(sf: SiteFibers, out_prefix: str, title: str = "", sort: str = "shuffle",
                   max_rows: int = 3000, seed: int = 0, smooth_bp: int = 11,
                   mnase: Optional[pd.DataFrame] = None, split_by_label: bool = True,
                   xlabel: str = "distance to site (bp)", dot_s: Optional[float] = None,
                   width: float = 6.0, formats=("png", "pdf"), dpi: int = 300,
                   threshold_note: str = ""):
    plt = style.setup_mpl()
    from matplotlib.gridspec import GridSpec
    h = sf.half
    labels = sf.labels if split_by_label else ["all"]
    groups = {}
    for lab in labels:
        idx = np.arange(len(sf.rows)) if lab == "all" else np.nonzero(sf.rows.label.values == lab)[0]
        if idx.size > max_rows:
            idx = np.sort(np.random.default_rng(seed).choice(idx, max_rows, replace=False))
        groups[lab] = order_rows(sf, idx, sort, seed)
    heights, kinds = [], []
    if mnase is not None and len(mnase):
        kinds.append("mnase"); heights.append(1.2)
    for lab in labels:
        kinds.append(("fib", lab)); heights.append(max(0.8, 3.5 * len(groups[lab]) / max(1, max_rows)))
    kinds.append("profile"); heights.append(1.4)
    fig = plt.figure(figsize=(width, sum(heights) + 0.8))
    gs = GridSpec(len(kinds), 1, figure=fig, height_ratios=heights, hspace=0.12)
    ax0 = None
    for r, k in enumerate(kinds):
        ax = fig.add_subplot(gs[r], sharex=ax0) if ax0 is not None else fig.add_subplot(gs[r])
        ax0 = ax0 or ax
        if k == "mnase":
            style.plot_density_scatter(mnase.rel.values, mnase.length.values, ax=ax, s=3)
            ax.axhline(147, color="gray", lw=0.6, ls=":")
            ax.set_ylim(0, mnase.length.max() + 5)
            ax.set_ylabel("MNase\nfrag (bp)", fontsize=style.LAB_FS)
        elif k == "profile":
            for i, lab in enumerate(sf.labels):
                n_rows = int((sf.rows.label == lab).sum())
                ax.plot(sf.x, smooth(sf.profile(lab), smooth_bp), lw=1.0, color=style.label_color(i),
                        label=f"{lab} ({n_rows:,} fibers)")
            ax.set_ylabel(f"m6A / called A\n({smooth_bp} bp mean)", fontsize=style.LAB_FS)
            ax.set_xlabel(xlabel, fontsize=style.LAB_FS)
            ax.legend(fontsize=style.LEG_FS, frameon=False)
            ax.axvline(0, color="k", lw=0.5, ls="--")
        else:
            lab = k[1]
            order = groups[lab]
            n = len(order)
            xs, ys = [], []
            for y, i in enumerate(order):
                m = sf.meth(i)
                xs.append(m); ys.append(np.full(m.size, y))
            s = dot_s if dot_s is not None else max(0.15, min(2.0, 300.0 / max(n, 1)))
            if xs:
                ax.scatter(np.concatenate(xs), np.concatenate(ys), s=s, c=style.METH_COLOR_ALT,
                           marker=".", edgecolors="none", rasterized=True)
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_ylabel(f"{lab}\n{n:,} fibers", fontsize=style.LAB_FS)
            ax.set_yticks([])
            ax.axvline(0, color="k", lw=0.4, ls="--", alpha=0.6)
        ax.tick_params(labelsize=style.TICK_FS, labelbottom=(k == "profile"))
    ax0.set_xlim(-h - 0.5, h + 0.5)
    if title:
        fig.suptitle(f"{title}  (sites used: {sf.n_sites_used:,}; sort: {sort})", fontsize=style.TITLE_FS)
    if threshold_note:
        fig.text(0.99, 0.003, f"m6A call: {threshold_note}", ha="right", fontsize=5, color="grey")
    paths = style.save(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return paths
