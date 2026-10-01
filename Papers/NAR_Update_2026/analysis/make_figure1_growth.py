#!/usr/bin/env python
"""Figure 1: what the database gained, and from where.

Regenerate with:  python make_figure1_growth.py
Shared palette, NUMBERS and helpers live in figure_common.py.
"""
from figure_common import *  # noqa: F401,F403 -- palette, NUMBERS, helpers
from figure_common import save


# ============================ FIGURE 1 ======================================
def figure1():
    """One row, three panels: what the database gained, and from where.

    A is the only panel on a count axis -- it compares 2020 with 2026, and a
    percentage axis would erase the thing it exists to show (KEGG flat at 17.8k
    while ChEBI arrives at 11.4k). B and C are 0-100% so the three sources are
    directly comparable; every count annotation was dropped in favour of the
    percentage, and the segment keys are direct-labelled on the top bar rather
    than carried in a legend outside the frame.
    """
    fig = plt.figure(figsize=(7.0, 1.95))
    gs = fig.add_gridspec(1, 3, left=0.107, right=0.984, top=0.985, bottom=0.135,
                          wspace=0.30)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))

    def key(ax, items):
        """Colour key inside the frame, upper right. Uses the legend layout
        engine rather than hand-placed swatches: two attempts at estimating
        text width in axes fractions put every swatch on top of its own label,
        because character width at 6pt is roughly twice what it looks like.
        Inside the frame because a legend under the axes is whitespace the row
        cannot afford."""
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(facecolor=c, edgecolor="none", label=n)
                           for n, c in items],
                  loc="upper right", frameon=False,
                  fontsize=6.0, ncol=len(items), handlelength=1.0, handleheight=0.85,
                  handletextpad=0.4, columnspacing=0.85, borderpad=0.1,
                  borderaxespad=0.98, labelcolor=INK2)

    def tag(ax, letter):
        ax.annotate(letter, xy=(0.012, 0.955), xycoords="axes fraction",
                    fontsize=9.5, fontweight="bold", va="top", ha="left", color=INK)

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
    top = len(rows) - 1
    a.text(400, top + h / 2 + 0.02, "2020", va="center", ha="left", fontsize=6.0,
           color=INK, zorder=5)
    a.text(400, top - h / 2 - 0.02, "2026", va="center", ha="left", fontsize=6.0,
           color=SURFACE, fontweight="bold", zorder=5)
    a.set_yticks([]); a.set_xlim(0, 30500); a.set_ylim(-0.6, len(rows) - 0.05)
    a.set_xticks([0, 10000, 20000, 30000]); a.set_xticklabels(["0", "10k", "20k", "30k"])
    strip(a); tag(a, "A")

    # -- B: share of each source's reactions that no other primary supplies.
    _by = {r[0]: r for r in NUMBERS["rxn_sources"]}
    rows = [_by[s] for s in ("MetaCyc", "KEGG", "Rhea")]
    ys = range(len(rows))[::-1]
    for i, (lab, total, uniq, _uc) in zip(ys, rows):
        pct = 100 * uniq / total if total else 0
        b.barh(i, pct, height=0.62, color=BLUE, zorder=4)
        b.barh(i, 100 - pct, left=pct, height=0.62, color=BLUE_200, zorder=3)
        b.text(pct - 1.6, i, f"{pct:.0f}%", va="center", ha="right", fontsize=6.6,
               color=SURFACE, fontweight="bold", zorder=5)
    key(b, [("unique", BLUE), ("shared", BLUE_200)])
    b.set_yticks([]); b.set_xlim(0, 100); b.set_ylim(-0.6, len(rows) - 0.05)
    b.set_xticks([0, 50, 100]); b.set_xticklabels(["0", "50", "100%"])
    strip(b); tag(b, "B")

    # -- C: of that unique contribution, the share whose every compound carries
    # a structure, so the reaction can be balanced and decomposed. Categories
    # are B's, in B's order, so they are labelled once.
    for i, (lab, _t, uniq, ucomp) in zip(ys, rows):
        pct = 100 * ucomp / uniq if uniq else 0
        c.barh(i, pct, height=0.62, color=BLUE, zorder=4)
        c.barh(i, 100 - pct, left=pct, height=0.62, color=NEUTRAL, zorder=3)
        c.text(pct - 1.6, i, f"{pct:.0f}%", va="center", ha="right", fontsize=6.6,
               color=SURFACE, fontweight="bold", zorder=5)
    key(c, [("complete", BLUE), ("incomplete", NEUTRAL)])
    c.set_yticks([]); c.set_xlim(0, 100); c.set_ylim(-0.6, len(rows) - 0.05)
    c.set_xticks([0, 50, 100]); c.set_xticklabels(["0", "50", "100%"])
    strip(c); tag(c, "C")
    return fig


if __name__ == "__main__":
    save(figure1(), "figure1_growth.pdf")
