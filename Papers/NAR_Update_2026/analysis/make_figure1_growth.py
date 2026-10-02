#!/usr/bin/env python
"""Figure 1: what the database gained, and from where.

Regenerate with:  python make_figure1_growth.py
Shared palette, NUMBERS and helpers live in figure_common.py. The pathway-class
panel that was D here is Supplementary Figure S1 (make_figureS1_pathway_classes.py).
"""
from figure_common import *  # noqa: F401,F403 -- palette, NUMBERS, helpers
from figure_common import save, _euler_regions


# ============================ FIGURE 1 ======================================
def figure1():
    """Three panels in one row: what the database gained, and from where.

    A is the only panel on a count axis -- it compares 2020 with 2026, and a
    percentage axis would erase the thing it exists to show (KEGG flat at 17.8k
    while ChEBI arrives at 11.4k). C is 0-100%. B is an Euler diagram rather
    than the stacked unique/shared bar it replaces: that bar reported what
    fraction of each source is unique but never said who the rest is shared
    WITH, which is the question integrating Rhea raises.

    The figure is drawn at the width it prints -- 7.0 in, the 178 mm text
    width of the OUP class, included at width=\\textwidth -- so every size here
    is the size on the page, and nothing is set below 7 pt. Every count is
    read from the released files through figure_common; none is transcribed.
    """
    fig = plt.figure(figsize=(7.0, 2.55))
    gs = fig.add_gridspec(1, 3, wspace=0.34, left=0.100, right=0.985,
                          top=0.955, bottom=0.150)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))
    for ax in (a, c):
        ax.tick_params(labelsize=7.4)

    def key(ax, items, **kw):
        """Colour key inside the frame. Uses the legend layout engine rather
        than hand-placed swatches: two attempts at estimating text width in
        axes fractions put every swatch on top of its own label, because
        character width at small sizes is roughly twice what it looks like."""
        from matplotlib.patches import Patch
        opts = dict(loc="upper right", frameon=False, fontsize=7.2,
                    ncol=len(items), handlelength=1.0, handleheight=0.85,
                    handletextpad=0.4, columnspacing=0.85, borderpad=0.1,
                    borderaxespad=0.98, labelcolor=INK2)
        opts.update(kw)
        ax.legend(handles=[Patch(facecolor=col, edgecolor="none", label=n)
                           for n, col in items], **opts)

    def tag(ax, letter, xy=(0.012, 0.955), va="top"):
        ax.annotate(letter, xy=xy, xycoords="axes fraction", fontsize=10.5,
                    fontweight="bold", va=va, ha="left", color=INK)

    # -- A: compounds per structure source, 2020 vs 2026. Percentage change only.
    rows = NUMBERS["growth"]; ys = range(len(rows))[::-1]; h = 0.40
    for i, (lab, old_, new_) in zip(ys, rows):
        if old_:
            a.barh(i + h / 2 + 0.02, old_, height=h, color=BLUE_200, zorder=3)
            change = 100 * (new_ - old_) / old_
            # "+-0%" is what a signed prefix on a rounded negative gives; a
            # change that rounds to nothing is written as such.
            pct = "±0%" if round(change) == 0 else f"{change:+.0f}%"
        else:
            pct = "new"
        a.barh(i - h / 2 - 0.02, new_, height=h, color=BLUE, zorder=3)
        a.text(-600, i, "ChEBI/Rhea" if lab == "ChEBI" else lab,
               va="center", ha="right", fontsize=7.8, color=INK)
        a.text(new_ + 500, i, pct, va="center", ha="left", fontsize=7.4,
               color=INK2, fontweight="bold")
    topi = len(rows) - 1
    a.text(450, topi + h / 2 + 0.02, "2020", va="center", ha="left", fontsize=7.0,
           color=INK, zorder=5)
    a.text(450, topi - h / 2 - 0.02, "2026", va="center", ha="left", fontsize=7.0,
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
    fig.text(bx.x0, bx.y1 - 0.004, "B", fontsize=10.5, fontweight="bold",
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
        c.text(pct - 1.8, i, f"{pct:.0f}%", va="center", ha="right", fontsize=7.4,
               color=SURFACE, fontweight="bold", zorder=5)
        c.text(-2.8, i, lab, va="center", ha="right", fontsize=7.8, color=INK)
    key(c, [("complete", BLUE), ("incomplete", NEUTRAL)])
    c.set_yticks([]); c.set_xlim(0, 100); c.set_ylim(-0.6, len(rows) - 0.05)
    c.set_xticks([0, 50, 100]); c.set_xticklabels(["0", "50", "100%"])
    strip(c); tag(c, "C")
    return fig


def _euler_panel(ax):
    """Three overlapping circles, one per primary reaction source.

    Fills are translucent so an overlap reads as the mix of its parents, with
    an opaque outline per set so membership stays legible in greyscale and
    under CVD -- identity is never carried by fill alone. The panel is about
    how the three sources overlap; reactions that come from none of them are
    not drawn.
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
        ax.text(x, y, f"{r[k]:,}", ha="center", va="center", fontsize=7.3,
                color=INK, zorder=5,
                fontweight="bold" if k in ("M", "K", "H") else "normal")
    # Set names ride outside their own circle, in the set's colour.
    for k, (x, y, ha) in {"M": (0.295, 0.978, "center"),
                          "K": (0.855, 0.978, "center"),
                          "H": (0.475, 0.018, "center")}.items():
        ax.text(x, y, f"{full[k]}  {tot[full[k]]:,}", ha=ha, va="center",
                fontsize=7.6, color=col[k], fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)


if __name__ == "__main__":
    save(figure1(), "figure1_growth.pdf")
