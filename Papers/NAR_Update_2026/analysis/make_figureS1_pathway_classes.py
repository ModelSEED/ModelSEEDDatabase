#!/usr/bin/env python
"""Supplementary Figure S1: where the growth since 2020 landed, by MetaCyc class.

Regenerate with:  python make_figureS1_pathway_classes.py
This was panel D of Figure 1 and moved to the supplement at the corresponding
author's request. It answers reviewer 1's question about where new biochemistry
is still arriving; shared palette and loaders live in figure_common.py.
"""
from figure_common import *  # noqa: F401,F403 -- palette, helpers
from figure_common import save, _pathway_classes

# Display names for the MetaCyc classes. The ontology ids are not reader-facing
# ("POLYKETIDE-SYN", "ALKALOIDS-SYN") and the full names are too long for an
# axis, so each is shortened once, here, rather than in the data function --
# which stays a faithful read of the source table.
CLASS_LABEL = {
    "Antibiotic-Biosynthesis": "antibiotics",
    "O-Antigen-Biosynthesis": "O-antigen",
    "POLYKETIDE-SYN": "polyketides",
    "Toxin-Biosynthesis": "toxins",
    "Lipid-Biosynthesis": "lipids",
    "Branched-Fatty-Acids-Biosynthesis": "branched fatty acids",
    "Fatty-acid-biosynthesis": "fatty acids",
    "ALKALOIDS-SYN": "alkaloids",
    "Sterol-Biosynthesis": "sterols",
    "Lipid-IV-A-Biosynthesis": "lipid IV-A",
    "ACYLSUGAR-BIOSYNTHESIS": "acylsugars",
    "Metabolic-Clusters": "metabolic clusters",
}


def _short(name):
    """Fall back to a trimmed ontology name when a class is not in the table."""
    return CLASS_LABEL.get(name, name.replace("-", " ").lower())


# ====================== SUPPLEMENTARY FIGURE S1 =============================
def figureS1():
    """Reactions gained since 2020 against those already held, by MetaCyc class.

    Two steps of one hue, because this is the same quantity at two times, not
    two identities -- the encoding Figure 1A uses for the compound sources.
    Drawn at 5.5 in for a 0.85\\textwidth include in the supplement, so the
    type is at or above 8 pt on the page.
    """
    fig = plt.figure(figsize=(5.5, 3.0))
    # Left margin is wide (room for "branched fatty acids" etc. at the bar
    # labels); right margin is mirrored to match it exactly, so the axes box
    # itself -- not just the figure canvas -- is centered on the page once
    # \includegraphics centers the whole PDF in S05's \figure environment.
    LEFT = 0.235
    gs = fig.add_gridspec(1, 1, left=LEFT, right=1 - LEFT, top=0.965, bottom=0.170)
    ax = fig.add_subplot(gs[0, 0])
    rows = _pathway_classes(top=8)
    ys = range(len(rows))[::-1]
    h = 0.40
    xmax = max(max(n, o) for _, n, o in rows) * 1.08
    for i, (name, new_, old_) in zip(ys, rows):
        ax.barh(i + h / 2 + 0.02, old_, height=h, color=BLUE_200, zorder=3)
        ax.barh(i - h / 2 - 0.02, new_, height=h, color=BLUE, zorder=3)
        ax.text(-xmax * 0.012, i, _short(name), va="center", ha="right",
                fontsize=8.4, color=INK)
    ax.set_yticks([]); ax.set_ylim(-0.62, len(rows) - 0.38)
    ax.set_xlim(0, xmax)
    ax.tick_params(labelsize=8.0)
    ax.set_xlabel("reactions", fontsize=8.4, color=INK2, labelpad=2.0)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=BLUE_200, edgecolor="none", label="already held in 2020"),
                       Patch(facecolor=BLUE, edgecolor="none", label="added since 2020")],
              loc="lower right", frameon=False, fontsize=7.8, ncol=1,
              handlelength=1.0, handleheight=0.85, handletextpad=0.4,
              labelspacing=0.3, borderpad=0.1, borderaxespad=0.75,
              labelcolor=INK2)
    strip(ax)
    # Outward x-ticks, overriding the shared Grace style's inward default for
    # this figure only (per corresponding-author request, 2026-10) -- the
    # y-axis carries no ticks (set_yticks([]) above), so only x needs it.
    ax.tick_params(axis="x", which="both", direction="out")
    return fig


if __name__ == "__main__":
    save(figureS1(), "figureS1_pathway_classes.pdf")
