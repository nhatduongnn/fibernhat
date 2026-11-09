"""Shared colours / fonts (values from the paper figure modules fig1/fig2_panels)."""
from __future__ import annotations

import numpy as np

METH_COLOR = "steelblue"
METH_COLOR_ALT = "#1f78b4"      # figure 1 aggregate / panel G blue
BACKBONE_COLOR = "lightgray"
STRAND_COLOR = {"+": "#1f3a5f", "-": "#8b1a1a", ".": "#888888"}
LABEL_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#7f7f7f"]
LAB_FS, TICK_FS, TITLE_FS, LEG_FS = 8, 7, 9, 6.5


def label_color(i: int) -> str:
    return LABEL_COLORS[i % len(LABEL_COLORS)]


def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none",
                         "font.family": ["DejaVu Sans"], "axes.spines.top": False,
                         "axes.spines.right": False})
    return plt


def plot_density_scatter(x, y, bw=(5, 15), cmap="magma_r", vmin=None, vmax=None,
                         ax=None, s=2, alpha=1.0, zorder=1, cbar=False):
    """Scatter coloured by 2-D KDE density (vendored from DMS_seq plot_utils)."""
    import matplotlib.pyplot as plt
    import statsmodels.api as sm
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if ax is None:
        ax = plt.gca()
    if len(x) == 0:
        return None
    try:
        kde = sm.nonparametric.KDEMultivariate(data=[x, y], var_type="cc", bw=bw)
        z = kde.pdf([x, y])
    except ValueError:
        z = np.zeros(len(x))
    o = np.argsort(z)
    x, y, z = x[o], y[o], z[o]
    ax.scatter(x, y, color="none", edgecolor="#c0c0c0", s=s, zorder=zorder, rasterized=True)
    sp = ax.scatter(x, y, c=z, s=s + 1, zorder=zorder, edgecolors="none", cmap=cmap,
                    rasterized=True, vmin=vmin, vmax=vmax, alpha=alpha)
    if cbar:
        plt.colorbar(sp, ax=ax)
    return sp


def save(fig, out_prefix: str, formats=("png", "pdf"), dpi: int = 300):
    import os
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    paths = []
    for f in formats:
        p = f"{out_prefix}.{f}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        paths.append(p)
    return paths
