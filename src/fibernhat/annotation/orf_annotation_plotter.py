"""Vendored from DMS_seq/src/script/orf_annotation_plotter.py (lab code). Only change: font family fallback."""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib import collections as mc
import matplotlib.cm as cm
import matplotlib.patheffects as path_effects


class ORFAnnotationPlotter:

    def __init__(self, orfs, introns):

        self.orfs = orfs
        self.introns = introns

        self.show_spines = True
        self.show_minor_ticks = False
        self.plot_tick_labels = True
        self.span = None
        self.chrom = None
        self.inset = (0, 10.0)
        self.triangle_width = 500
        self.text_horizontal_offset = 10
        self.plot_orf_names = False
        self.height = 40
        self.text_vertical_offset = 18
        self.y_padding = 10

    def set_span_chrom(self, span, chrom):
        if not self.span == span or not chrom == self.chrom:
            self.chrom = chrom
            self.span = int(span[0]), int(span[1])

        span_width = self.span[1] - self.span[0]
        self.span_width = span_width

        # triangle is proportion of the span width
        self.triangle_width = span_width/100.

        # overlap triangle with rect
        self.epsilon = span_width * 1e-4

        self.text_horizontal_offset = span_width * 1e-2

    def get_color(self, watson, gene_type):
        colors = {
            "Verified": cm.Blues(0.8),
            "Putative": cm.Blues(0.5),
            "Dubious": "#898888"
        }
        if not watson:
            colors['Verified'] = cm.Reds(0.5)
            colors['Putative'] = cm.Reds(0.3)

        colors['Uncharacterized'] = colors['Putative']
        color = colors[gene_type]
        return color

    def plot_orf_annotation(self, ax, start, end, name, CDS_introns=None, 
        watson=True, offset=False, gene_type="Verified", TSS=np.nan, PAS=np.nan,
        plot_orf_names=True, text_horizontal_offset=None, y=None, color=None, 
        text_style='italic', fontsize=12):

        if color is None:
            color = self.get_color(watson, gene_type)

        rect_width = end - start - self.triangle_width

        if text_horizontal_offset is None:
            text_horizontal_offset = self.text_horizontal_offset

        inset = self.inset

        if watson:
            rect_start = start
            triangle_start = end - self.triangle_width - self.epsilon*2 # add one to overlap triangle with rect
            y_baseline = offset * (self.y_padding + self.height)
            text_start = start + text_horizontal_offset

            if np.isnan(TSS): TSS = start
            if np.isnan(PAS): PAS = end
        else:
            triangle_start = start + self.epsilon # add one to overlap triangle with rect
            rect_start = start + self.triangle_width
            y_baseline = (offset+1) * (-self.y_padding-self.height)
            text_start = end - text_horizontal_offset
            if np.isnan(TSS): TSS = end
            if np.isnan(PAS): PAS = start

        if y is not None:
            y_baseline = y

        # plot pointed rectangle for ORF
        plot_rect(ax, rect_start, y_baseline, rect_width, self.height, color, 
            inset=inset)

        plot_iso_triangle(ax, triangle_start, y_baseline, 
            self.triangle_width, self.height, color, facing_right=watson,
            inset=inset[1])

        # plot introns as rectangles overtop of ORF
        if CDS_introns is not None and len(CDS_introns) > 0:
            for idx, intron in CDS_introns.iterrows():
                if intron['cat'] == 'intron':
                    intron_width = intron.stop - intron.start
                    # plot introns as a lighter version of the CDS
                    plot_rect(ax, intron.start, y_baseline, intron_width, 
                        self.height, 'white', inset=inset)
                    plot_rect(ax, intron.start, y_baseline, intron_width, 
                        self.height, color, fill_alpha=0.25, inset=inset)

        # plot name of ORF
        if plot_orf_names:
            plot_text(ax,
                text_start, 
                y_baseline+self.text_vertical_offset, 
                name, color, flipped=(not watson),
                fontsize=fontsize,
                text_style=text_style)

        # Plot TSS PAS line indicators
        if TSS != np.nan or PAS != np.nan:
            plot_TSS_PAS(ax, start, end, TSS, PAS, 
                y_baseline, self.height, color, flipped=watson, inset=inset[1])

        draw_TSS_arrow(ax, TSS, y_baseline+inset[1], flip=not watson, 
            plot_width=self.triangle_width, color=color)

    def plot_orf_annotations(self, ax, 
        orf_classes=['Verified', 'Uncharacterized', 'Dubious'],
        custom_orfs=None, should_auto_offset=True):

        if custom_orfs is None:
            orfs = self.orfs
        else: orfs = custom_orfs

        span = self.span
        chrom = self.chrom
        genes = orfs[(orfs['chr'] == int(chrom)) & 
                          (orfs['stop'] > span[0]) & 
                          (orfs['start'] < span[1]) & 
                          (orfs.orf_class.isin(orf_classes))]

        try:
            genes = genes.sort_values(['strand', 'start']).reset_index()

        # older pandas version
        except AttributeError:
            genes = genes.sort(['strand', 'start']).reset_index()

        ax.set_ylim(-100, 100)

        tick_intervals = 500, 100
        
        for idx, gene in genes.iterrows():

            name = gene['orf_name']

            if self.plot_orf_names:
                if not name == gene['name']:
                    name = "{}/{}".format(gene['name'], gene['orf_name'])
            else:
                name = gene['name']

            if not self.introns is None:
                gene_introns = self.introns[
                    self.introns['parent'] == gene['orf_name']]
            else: gene_introns = None

            offset = False
            if should_auto_offset: offset = (idx % 2) == 0

            self.plot_orf_annotation(ax, gene.start,
                      gene.stop, name,
                      CDS_introns=gene_introns,
                      watson=(gene.strand == "+"),
                      gene_type=gene.orf_class,
                      offset=offset,
                      TSS=gene['TSS'],
                      PAS=gene['PAS'])


        ax.set_xticks(np.arange(span[0], span[1]+500, 2000), minor=False)
        ax.set_xticks(np.arange(span[0], span[1]+500, 1000), minor=True)
        ax.set_xlim(span[0], span[1])
        ax.set_yticks([])

        return ax

    def plot_centromere(self, ax, start, end, name,
            watson=True, offset=False,
            text_horizontal_offset=None, y=None):

            color='#4A4A4A'

            rect_width = end - start - self.triangle_width

            if text_horizontal_offset is None:
                text_horizontal_offset = self.text_horizontal_offset

            inset = self.inset

            if watson:
                rect_start = start
                triangle_start = end - self.triangle_width - self.epsilon*2 # add one to overlap triangle with rect
                y_baseline = offset * (self.y_padding + self.height)
                text_start = start
                ha = 'left'
            else:
                triangle_start = start + self.epsilon # add one to overlap triangle with rect
                rect_start = start + self.triangle_width
                y_baseline = (offset+1) * (-self.y_padding-self.height)
                text_start = end
                ha = 'right'

            if y is not None:
                y_baseline = y

            # plot pointed rectangle
            plot_rect(ax, rect_start, y_baseline, rect_width, self.height, color, 
                inset=inset)

            plot_iso_triangle(ax, triangle_start, y_baseline, 
                self.triangle_width, self.height, color, facing_right=watson,
                inset=inset[1])

            text = ax.text(text_start, y_baseline+self.text_vertical_offset,
                name, fontsize=12, clip_on=True, zorder=65, 
                           rotation=0, va='center', ha=ha,
                           fontdict={'family': ['Open Sans', 'DejaVu Sans']},
                           color='white')
            text.set_path_effects([path_effects.Stroke(linewidth=1.5, foreground=color),
                                   path_effects.Normal()])


def plot_TSS_PAS(ax, start, end, TSS, PAS, y, height, color, flipped=False, 
    inset=0):
    """# Plot TSS PAS lines, if TSS or PAS does not exist use ORF start
        # and end boundaries"""

    y_bottom = y+inset/2.
    y_top = y+height-inset/2

    if TSS is not None:
        ax.plot([TSS, TSS], [y_bottom, y_top], color=color, lw=1.5, zorder=111)

    if PAS is not None:
        ax.plot([PAS, PAS], [y_bottom, y_top], color=color, lw=1.5, zorder=111)

    if not flipped:
        if np.isnan(TSS): TSS = start
        if np.isnan(PAS): PAS = end
    else:
        if np.isnan(TSS): TSS = end
        if np.isnan(PAS): PAS = start

    TSS_span = TSS, PAS    
    ax.plot(TSS_span, [y+height/2., y+height/2.], lw=1.5, color=color, zorder=60)


def draw_TSS_arrow(ax, x, y, flip, plot_width, color):

    tri_width, tri_height = 45, 10
    tss_height = 33.
    tss_width = 12
    lw=1.5

    if flip: 
        tss_width = -tss_width
        tri_width = -tri_width

    plot_iso_triangle(ax, x+tss_width, y+tss_height-tri_height/2, tri_width, 
        tri_height, color)

    ax.plot([x, x+tss_width], [y+tss_height, y+tss_height], c=color, lw=lw)
    ax.plot([x, x], [y, y+tss_height], c=color, lw=lw)


def plot_iso_triangle(ax, x, y, width, height, color, facing_right=True,
    inset=0):
    """
    Plot left or right pointing isosceles triangle for end cap of ORF
    """

    y_bottom = y+inset/2.
    y_top = y+height-inset/2.
    x_flat = x

    if facing_right: 
        x_pointed = x+width
    else: 
        x_pointed = x
        x_flat += width

    X = np.array([[x_flat,y_bottom], [x_flat,y_top], 
        [x_pointed, (y_bottom+y_top)/2.0]])

    triangle = plt.Polygon(X, color=color, linewidth=0, zorder=2)
    ax.add_patch(triangle)

def plot_text(ax, start, y, name, color, flipped=False, text_style='italic',
              fontsize=22, rotation='horizontal', ha=None):

    if ha is None:
        ha = 'left'
        if flipped: ha = 'right'

    text = ax.text(start, y, name, fontsize=fontsize, clip_on=True, zorder=65, 
                   rotation=rotation, va='center', ha=ha, 
                   fontdict={'family': ['Open Sans', 'DejaVu Sans'], 'fontweight': 'regular',
                   'style': text_style},
                   color='white')
    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground=color),
                           path_effects.Normal()])


def plot_rect(ax, x, y, width, height, color=None, facecolor=None, 
    edgecolor=None, ls='solid', fill_alpha=1., zorder=40, lw=0.0, 
    inset=(0.0, 0.0), fill=True, joinstyle='round'):
    """
    Plot a rectangle for ORF plotting
    """

    if edgecolor is None: edgecolor = color
    if facecolor is None: facecolor = color

    patch = ax.add_patch(
                    patches.Rectangle(
                        (x, y + inset[1]/2.0),   # (x,y)
                        width, height - inset[1], # size
                        facecolor=facecolor,
                        edgecolor=edgecolor,
                        lw=lw,
                        joinstyle=joinstyle,
                        ls=ls,
                        fill=fill,
                        alpha=fill_alpha,
                        zorder=zorder
                    ))
