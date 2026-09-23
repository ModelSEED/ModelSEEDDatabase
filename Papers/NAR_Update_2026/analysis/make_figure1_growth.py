#!/usr/bin/env python
"""Figure 1: what the database gained, and from where.

Regenerate with:  python make_figure1_growth.py
Shared palette, NUMBERS and helpers live in figure_common.py.
"""
from figure_common import *  # noqa: F401,F403 -- palette, NUMBERS, helpers
from figure_common import save, _euler_regions, _pathway_classes

# Display names for the MetaCyc classes in panel D. The ontology ids are not
# reader-facing ("POLYKETIDE-SYN", "ALKALOIDS-SYN") and the full names are too
# long for an axis at 6pt, so each is shortened once, here, rather than in the
# data function -- which stays a faithful read of the source table.
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


# ============================ FIGURE 1 ======================================
def figure1():
    """Four panels: what the database gained, and from where.

    A is the only panel on a count axis -- it compares 2020 with 2026, and a
    percentage axis would erase the thing it exists to show (KEGG flat at 17.8k
    while ChEBI arrives at 11.4k). C is 0-100%. B is an Euler diagram rather
    than the stacked unique/shared bar it replaces: that bar reported what
    fraction of each source is unique but never said who the rest is shared
    WITH, which is the question integrating Rhea raises. D answers reviewer 1's
    question about where new biochemistry is still arriving.

    Every count is read from the released files through figure_common; none is
    transcribed.
    """
    fig = plt.figure(figsize=(7.0, 3.42))
    # Two gridspecs rather than one with a nested row: panel D's class labels
    # are longer than anything in the top row and need their own left inset.
    top = fig.add_gridspec(1, 3, wspace=0.30, left=0.105, right=0.984,
                           top=0.978, bottom=0.560)
    bot = fig.add_gridspec(1, 1, left=0.152, right=0.984, top=0.432, bottom=0.096)
    a, b, c = (fig.add_subplot(top[0, i]) for i in range(3))
    d = fig.add_subplot(bot[0, 0])

    def key(ax, items, **kw):
        """Colour key inside the frame. Uses the legend layout engine rather
        than hand-placed swatches: two attempts at estimating text width in
        axes fractions put every swatch on top of its own label, because
        character width at 6pt is roughly twice what it looks like."""
        from matplotlib.patches import Patch
        opts = dict(loc="upper right", frameon=False, fontsize=6.0,
                    ncol=len(items), handlelength=1.0, handleheight=0.85,
                    handletextpad=0.4, columnspacing=0.85, borderpad=0.1,
                    borderaxespad=0.98, labelcolor=INK2)
        opts.update(kw)
        ax.legend(handles=[Patch(facecolor=col, edgecolor="none", label=n)
                           for n, col in items], **opts)

    def tag(ax, letter, xy=(0.012, 0.955), va="top"):
        ax.annotate(letter, xy=xy, xycoords="axes fraction", fontsize=9.5,
                    fontweight="bold", va=va, ha="left", color=INK)

    # -- A: compounds per structure source, 2020 vs 2026. Percentage change only.
    rows = NUMBERS["growth"]; ys = range(len(rows))[::-1]; h = 0.40
    for i, (lab, old_, new_) in zip(ys, rows):
        if old_:
            a.barh(i + h / 2 + 0.02, old_, height=h, color=BLUE_200, zorder=3)
            pct = f"+{100*(new_-old_)/old_:.0f}%"
        else:
            pct = "new"
        a.barh(i - h / 2 - 0.02, new_, height=h, color=BLUE, zorder=3)
        a.text(-500, i, "ChEBI/Rhea" if lab == "ChEBI" else lab,
               va="center", ha="right", fontsize=7.0, color=INK)
        a.text(new_ + 450, i, pct, va="center", ha="left", fontsize=6.8,
               color=INK2, fontweight="bold")
    topi = len(rows) - 1
    a.text(400, topi + h / 2 + 0.02, "2020", va="center", ha="left", fontsize=6.0,
           color=INK, zorder=5)
    a.text(400, topi - h / 2 - 0.02, "2026", va="center", ha="left", fontsize=6.0,
           color=SURFACE, fontweight="bold", zorder=5)
    a.set_yticks([]); a.set_xlim(0, 30500); a.set_ylim(-0.6, len(rows) - 0.05)
    a.set_xticks([0, 10000, 20000, 30000]); a.set_xticklabels(["0", "10k", "20k", "30k"])
    strip(a); tag(a, "A")

    # -- B: three-set Euler over reactions. Circle areas are not solvable
    # exactly for three sets, so the geometry is the conventional symmetric
    # arrangement and every region carries its count -- the counts are the
    # quantitative channel, the circles only carry membership.
    _euler_panel(b)
    bx = b.get_subplotspec().get_position(fig)
    fig.text(bx.x0, bx.y1 - 0.004, "B", fontsize=9.5, fontweight="bold",
             va="top", ha="left", color=INK)

    # -- C: share of each source's UNIQUE contribution that is complete, meaning
    # every reagent carries a structure so the reaction can be balanced and
    # decomposed. Rows are B's sources, in B's order.
    _by = {r[0]: r for r in NUMBERS["rxn_sources"]}
    rows = [_by[s] for s in ("MetaCyc", "KEGG", "Rhea")]
    ys = range(len(rows))[::-1]
    for i, (lab, _t, uniq, ucomp) in zip(ys, rows):
        pct = 100 * ucomp / uniq if uniq else 0
        c.barh(i, pct, height=0.62, color=BLUE, zorder=4)
        c.barh(i, 100 - pct, left=pct, height=0.62, color=NEUTRAL, zorder=3)
        c.text(pct - 1.6, i, f"{pct:.0f}%", va="center", ha="right", fontsize=6.6,
               color=SURFACE, fontweight="bold", zorder=5)
        c.text(-2.5, i, lab, va="center", ha="right", fontsize=6.6, color=INK)
    key(c, [("complete", BLUE), ("incomplete", NEUTRAL)])
    c.set_yticks([]); c.set_xlim(0, 100); c.set_ylim(-0.6, len(rows) - 0.05)
    c.set_xticks([0, 50, 100]); c.set_xticklabels(["0", "50", "100%"])
    strip(c); tag(c, "C")

    # -- D: where the growth landed, by MetaCyc class.
    _pathway_panel(d, key)
    tag(d, "D", xy=(0.0, 1.015), va="bottom")
    return fig


def _euler_panel(ax):
    """Three overlapping circles, one per primary reaction source.

    Fills are translucent so an overlap reads as the mix of its parents, with
    an opaque outline per set so membership stays legible in greyscale and
    under CVD -- identity is never carried by fill alone.
    """
    from matplotlib.patches import Circle
    r = _euler_regions()
    tot = r["totals"]
    full = {"M": "MetaCyc", "K": "KEGG", "H": "Rhea"}
    # Conventional symmetric three-circle layout.
    R = 0.335
    cen = {"M": (0.375, 0.600), "K": (0.625, 0.600), "H": (0.500, 0.390)}
    col = {"M": BLUE, "K": ORANGE, "H": AQUA}
    for k in ("M", "K", "H"):
        ax.add_patch(Circle(cen[k], R, facecolor=col[k], alpha=0.28,
                            edgecolor="none", zorder=2))
        ax.add_patch(Circle(cen[k], R, facecolor="none", edgecolor=col[k],
                            linewidth=0.9, zorder=4))
    # Region counts, at the visual centroid of each lens.
    at = {"M": (0.200, 0.690), "K": (0.800, 0.690), "H": (0.500, 0.170),
          "MK": (0.500, 0.755), "MH": (0.320, 0.395), "KH": (0.680, 0.395),
          "MKH": (0.500, 0.530)}
    for k, (x, y) in at.items():
        ax.text(x, y, f"{r[k]:,}", ha="center", va="center", fontsize=5.9,
                color=INK, zorder=5,
                fontweight="bold" if k in ("M", "K", "H") else "normal")
    # Set names ride outside their own circle, in the set's colour.
    for k, (x, y, ha) in {"M": (0.295, 0.978, "center"),
                          "K": (0.855, 0.978, "center"),
                          "H": (0.475, 0.018, "center")}.items():
        ax.text(x, y, f"{full[k]}  {tot[full[k]]:,}", ha=ha, va="center",
                fontsize=6.2, color=col[k], fontweight="bold")
    ax.text(1.045, 0.018, f"+{r['none']:,} from other sources", ha="right",
            va="center", fontsize=5.5, color=MUTED)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)


def _pathway_panel(ax, key):
    """Reactions gained since 2020 against those already held, by MetaCyc class.

    Two steps of one hue, because this is the same quantity at two times, not
    two identities -- the same encoding panel A uses for the compound sources.
    """
    rows = _pathway_classes(top=8)
    ys = range(len(rows))[::-1]
    h = 0.40
    xmax = max(max(n, o) for _, n, o in rows) * 1.22
    for i, (name, new_, old_) in zip(ys, rows):
        ax.barh(i + h / 2 + 0.02, old_, height=h, color=BLUE_200, zorder=3)
        ax.barh(i - h / 2 - 0.02, new_, height=h, color=BLUE, zorder=3)
        ax.text(-xmax * 0.012, i, _short(name), va="center", ha="right",
                fontsize=6.2, color=INK)
        ax.text(max(new_, old_) + xmax * 0.012, i, f"{new_:,} new", va="center",
                ha="left", fontsize=6.0, color=INK2, fontweight="bold")
    ax.set_yticks([]); ax.set_ylim(-0.62, len(rows) - 0.38)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("reactions", fontsize=6.4, color=INK2, labelpad=1.5)
    key(ax, [("already held in 2020", BLUE_200), ("added since 2020", BLUE)],
        loc="lower right", borderaxespad=0.75)
    strip(ax)


if __name__ == "__main__":
    save(figure1(), "figure1_growth.pdf")
