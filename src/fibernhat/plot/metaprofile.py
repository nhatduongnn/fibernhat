"""Q1..Q5 metaprofiles (viridis), one row per label, one column per site class."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import style

TITLES = {"dyad": "Nucleosome dyads", "nuc": "Nucleosome dyads", "tss": "TSS",
          "abf1": "Abf1 motifs", "reb1": "Reb1 motifs"}
XLAB = {"dyad": "distance to dyad (bp)", "nuc": "distance to dyad (bp)", "tss": "distance to TSS (bp)"}


def plot_quintile_profiles(prof: pd.DataFrame, summ: pd.DataFrame, out_prefix: str,
                           smooth_bp: int = 21, normalize: bool = True,
                           formats=("png", "pdf"), dpi: int = 200, threshold_note: str = "", show_note: bool = False):
    plt = style.setup_mpl()
    labels = list(dict.fromkeys(summ.label))
    classes = list(dict.fromkeys(prof.site_class))
    fig, axes = plt.subplots(len(labels), len(classes), figsize=(4.2 * len(classes), 2.8 * len(labels)),
                             sharex="col", squeeze=False)
    cmap = plt.get_cmap("viridis")
    for r, lab in enumerate(labels):
        for c, cls in enumerate(classes):
            ax = axes[r][c]
            for q in range(1, 6):
                d = prof[(prof.label == lab) & (prof.quintile == q) & (prof.site_class == cls)]
                if d.empty:
                    continue
                g = summ[(summ.label == lab) & (summ.quintile == q)].genome_m6a_rate.iloc[0]
                y = d.m6a_rate.values / g if normalize else d.m6a_rate.values
                y = pd.Series(y).rolling(smooth_bp, center=True, min_periods=1).mean().values
                lbl = f"Q{q}" + (" (least m6A)" if q == 1 else " (most)" if q == 5 else "")
                ax.plot(d.pos, y, color=cmap((q - 1) / 4), lw=1.2, label=lbl)
            if normalize:
                ax.axhline(1, color="#8a8f96", lw=0.6)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)   # as plot_quintile_pair.py L39
            if r == 0:
                n = prof[(prof.site_class == cls)].n_sites.max()
                ax.set_title(f"{TITLES.get(cls, cls)} (n={n:,})", fontsize=10)
            if r == len(labels) - 1:
                ax.set_xlabel(XLAB.get(cls, "distance to site (bp)"), fontsize=9)
            if c == 0:
                ax.set_ylabel(f"{lab}\n" + ("m6A / quintile genome mean" if normalize else "m6A fraction"),
                              fontsize=8)
    axes[0][-1].legend(fontsize=7, frameon=False, title="per-read m6A quintile", title_fontsize=7)
    if threshold_note and show_note:
        fig.text(0.99, 0.002, f"m6A call: {threshold_note}", ha="right", fontsize=6, color="grey")
    fig.tight_layout()
    paths = style.save(fig, out_prefix, formats, dpi)
    plt.close(fig)
    return paths
