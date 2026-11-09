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


AGG_PANELS = ("mnase", "fibers", "counts", "ratio")
AGG_HEIGHTS = {"mnase": 1.2, "counts": 1.0, "ratio": 1.0}


def aggregate_plot(sf: SiteFibers, out_prefix: str, title: str = "", sort: str = "meth_in_win",
                   max_rows: int = 3000, seed: int = 0, smooth_bp: int = 1,
                   mnase: Optional[pd.DataFrame] = None, split_by_label: bool = True,
                   xlabel: str = "distance to site (bp)", dot_s: Optional[float] = None,
                   panels=AGG_PANELS, width: float = 6.0, formats=("png", "pdf"), dpi: int = 300,
                   threshold_note: str = "", show_note: bool = False):
    """panels, top to bottom, any subset/order of: mnase, fibers (one panel per
    sample unless split_by_label=False), counts ("Methylation pileup"),
    ratio ("Methylated As / Possible As"). Default follows
    Fiber_seq_plots.process_aggregate_fiber_plot (fibers, counts, ratio; ratio
    unsmoothed) with the MNase V-plot on top when available."""
    plt = style.setup_mpl()
    from matplotlib.gridspec import GridSpec
    panels = [x.strip() for x in (panels.split(",") if isinstance(panels, str) else panels) if x.strip()]
    bad = [x for x in panels if x not in AGG_PANELS]
    if bad:
        raise ValueError(f"unknown panel(s) {bad}; choose from {','.join(AGG_PANELS)}")
    h = sf.half
    labels = sf.labels if split_by_label else ["all"]
    groups = {}
    for lab in labels:
        idx = np.arange(len(sf.rows)) if lab == "all" else np.nonzero(sf.rows.label.values == lab)[0]
        if idx.size > max_rows:
            idx = np.sort(np.random.default_rng(seed).choice(idx, max_rows, replace=False))
        groups[lab] = order_rows(sf, idx, sort, seed)
    heights, kinds = [], []
    for pn in panels:
        if pn == "mnase" and (mnase is None or not len(mnase)):
            continue
        if pn == "fibers":
            for lab in labels:
                kinds.append(("fib", lab)); heights.append(max(1.1, 3.5 * len(groups[lab]) / max(1, max_rows)))
        else:
            kinds.append(pn); heights.append(AGG_HEIGHTS[pn])
    fig = plt.figure(figsize=(width, sum(heights) + 0.8))
    gs = GridSpec(len(kinds), 1, figure=fig, height_ratios=heights, hspace=0.12)
    ax0 = None
    many = len(sf.labels) > 1
    for r, k in enumerate(kinds):
        ax = fig.add_subplot(gs[r], sharex=ax0) if ax0 is not None else fig.add_subplot(gs[r])
        ax0 = ax0 or ax
        if k == "mnase":
            style.plot_density_scatter(mnase.rel.values, mnase.length.values, ax=ax, s=3)
            ax.axhline(147, color="gray", lw=0.6, ls=":")
            ax.set_ylim(0, mnase.length.max() + 5)
            ax.set_ylabel("MNase fragment\nlength (bp)", fontsize=style.LAB_FS)
        elif k in ("counts", "ratio"):
            for i, lab in enumerate(sf.labels):
                n_rows = int((sf.rows.label == lab).sum())
                y = sf.mod[lab] if k == "counts" else sf.profile(lab)
                ax.plot(sf.x, smooth(y, smooth_bp), lw=1.0, color=style.label_color(i),
                        label=f"{lab} ({n_rows:,} fibers)")
            if k == "counts":
                ax.set_ylabel("Methylation\npileup", fontsize=style.LAB_FS)
            else:
                ax.set_ylabel("Methylated As /\nPossible As" + (f"\n({smooth_bp} bp mean)" if smooth_bp > 1 else ""),
                              fontsize=style.LAB_FS)
            if many:
                ax.legend(fontsize=style.LEG_FS, frameon=False)
        else:
            lab = k[1]
            order = groups[lab]
            n = len(order)
            xs, ys = [], []
            for y, i in enumerate(order):
                m = sf.meth(i)
                xs.append(m); ys.append(np.full(m.size, y))
            s_ = dot_s if dot_s is not None else max(0.15, min(2.0, 300.0 / max(n, 1)))
            if xs:
                ax.scatter(np.concatenate(xs), np.concatenate(ys), s=s_, c=style.METH_COLOR_ALT,
                           marker=".", edgecolors="none", rasterized=True)
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_ylabel(f"{lab} fibers" if many else "Number of fibers", fontsize=style.LAB_FS)
        ax.tick_params(labelsize=style.TICK_FS, labelbottom=(r == len(kinds) - 1))
        if r == 0 and title:
            ax.set_title(title, fontsize=style.TITLE_FS)
    fig.axes[-1].set_xlabel(xlabel, fontsize=style.LAB_FS)
    ax0.set_xlim(-h - 0.5, h + 0.5)
    if threshold_note and show_note:
        fig.text(0.99, 0.003, f"m6A call: {threshold_note}", ha="right", fontsize=5, color="grey")
    paths = style.save(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return paths
