#!/usr/bin/env python
"""Figure 3: direction assigned per source, agreement, and what underpins each grade.

Regenerate with:  python make_figure3_direction.py
Shared palette, NUMBERS and helpers live in figure_common.py.
"""
from figure_common import *  # noqa: F401,F403 -- palette, NUMBERS, helpers
from figure_common import save


# ============================ FIGURE 3 ======================================
def figure3():
    """Three panels in one row (2026-09-15).

    Was A/B/C with a six-lane C spanning the bottom. B (eQuilibrator against
    dGPredictor) dropped at Sam's request; the six-lane panel split into two
    equal panels, one per axis, which also retires the two-legends-in-one-axes
    hack -- the palette is reused across the axes (BLUE is self-certain here and
    corroborated there), so a shared key was never safe. Separate panels give
    each its own.
    """
    fig = plt.figure(figsize=(7.0, 2.30))
    # B and C carry no y labels and A's are abbreviated, so the panels run to
    # the page edges; wspace is the only furniture left between them.
    gs = fig.add_gridspec(1, 3, left=0.052, right=0.998, top=0.955,
                          bottom=0.150, wspace=0.155)
    a = fig.add_subplot(gs[0, 0])
    c = fig.add_subplot(gs[0, 1]); d = fig.add_subplot(gs[0, 2])
    # Gaps are asymmetric and gridspec wspace is not, so place by hand: B needs
    # room on its left for the gold/silver/bronze labels it carries for both
    # itself and C; C needs none, so it sits tight against B.
    L, R, GAP_AB, GAP_BC = 0.052, 0.998, 0.058, 0.016
    W = (R - L - GAP_AB - GAP_BC) / 3.0
    for _ax, _x0 in ((a, L), (c, L + W + GAP_AB), (d, L + 2 * W + GAP_AB + GAP_BC)):
        _b = _ax.get_position()
        _ax.set_position([_x0, _b.y0, W, _b.height])

    def tag(ax, letter, x=0.017, y=0.985):
        ax.annotate(letter, xy=(x, y), xycoords="axes fraction", fontsize=9.0,
                    fontweight="bold", va="top", ha="left", color=INK, zorder=6)

    FWD, REV, BACK = "#1c5cab", NEUTRAL, "#c2410c"
    UND = NEUTRAL
    # abbreviated so the y labels cost almost no width; the caption expands them
    SHORT_SRC = {"eQuilibrator": "eQ", "Group contribution": "GC",
                 "dGPredictor": "dG", "LLMs": "LLMs"}
    rows = NUMBERS["direction"]; ys = range(len(rows))[::-1]
    AMAX = max(f + e + r + q for _, f, e, r, q in rows)
    for i, (lab, f, e, r, q) in zip(ys, rows):
        tot = f + e + r + q; x = 0
        for v, col, nm in ((f, FWD, "\u2192"), (e, AQUA, "\u2194"),
                           (r, BACK, "\u2190"), (q, UND, "?")):
            pct = v
            a.barh(i, pct, left=x, height=0.66, color=col, zorder=3,
                   edgecolor=FRAME, lw=0.45)
            x += pct     # no in-bar glyphs (2026-09-15); the key carries them
        a.text(-AMAX * 0.028, i, SHORT_SRC.get(lab, lab),
               va="center", ha="right", fontsize=7.2, color=INK)
    a.set_xlim(0, AMAX); a.set_ylim(-0.60, len(rows) - 0.02); a.set_yticks([])
    a.set_xticks([0, 20000, 40000]); a.set_xticklabels(["0", "20k", "40k"])
    strip(a)
    # key as a COLUMN in the white space right of the shortest row (eQuilibrator,
    # 25k against an axis running to 46k), not a four-across strip in the
    # headroom, which crowded the panel letter
    from matplotlib.patches import Patch as _Patch
    a.legend(handles=[_Patch(facecolor=cc, edgecolor="none", label=nn)
                      for nn, cc in (("forward", FWD), ("reversible", AQUA),
                                     ("reverse", BACK), ("undet.", UND))],
             loc="upper right", bbox_to_anchor=(1.0, 0.915), frameon=False,
             fontsize=5.8, ncol=1, handlelength=0.9, handleheight=0.8,
             handletextpad=0.32, labelspacing=0.30, borderpad=0.1,
             borderaxespad=0.35, labelcolor=INK2)
    tag(a, "A")

    # C and D -- the two axes a grade is built from, one panel each, equal size.
    GR = ["gold", "silver", "bronze"]
    ORDER = {"grade_assess": ["measured", "self-certain", "self-confident",
                              "unconfident"],
             "grade_cross":  ["corroborated", "disputed", "unpaired",
                              "neither way"]}
    DMAX = max(sum(v for _, v, _ in NUMBERS[k].get(g, [])) or 1
               for k in ("grade_assess", "grade_cross") for g in GR)
    for ax, key, letter, title in ((c, "grade_assess", "B", "by self-assessment"),
                                   (d, "grade_cross", "C", "by cross-source")):
        lanes = [(g, NUMBERS[key].get(g, [])) for g in GR]
        yy = range(len(lanes))[::-1]
        seen = []   # filled in encounter order, then sorted to Table 1's below
        for i, (lab, segs) in zip(yy, lanes):
            x = 0
            for nm, v, col in segs:
                ax.barh(i, v, left=x, height=0.62, color=col, zorder=3,
                        edgecolor=FRAME, lw=0.45)
                if nm not in [n for n, _ in seen]:
                    seen.append((nm, col))
                x += v
            if key == "grade_assess":          # C repeats B's lanes; label once
                ax.text(-DMAX * 0.030, i, lab, va="center", ha="right",
                        fontsize=7.0, color=INK)
        ax.set_xlim(0, DMAX); ax.set_ylim(-0.62, len(lanes) - 0.32)
        ax.set_yticks([])
        ax.set_xticks([0, 5000, 10000, 15000])
        ax.set_xticklabels(["0", "5k", "10k", "15k"])
        strip(ax)
        # same anchor as A's key, so all three sit on one line across the figure
        rank = {n_: i for i, n_ in enumerate(ORDER[key])}
        seen.sort(key=lambda t: rank.get(t[0], 99))
        ax.legend(handles=[_Patch(facecolor=c_, edgecolor="none", label=n_)
                           for n_, c_ in seen],
                  loc="upper right", bbox_to_anchor=(1.0, 0.915), frameon=False,
                  fontsize=5.8, ncol=1, handlelength=0.9, handleheight=0.8,
                  handletextpad=0.32, labelspacing=0.30, borderpad=0.1,
                  borderaxespad=0.35, labelcolor=INK2)
        ax.text(0.5, -0.20, title, transform=ax.transAxes, ha="center",
                va="top", fontsize=6.4, color=MUTED)
        tag(ax, letter)
    return fig


if __name__ == "__main__":
    save(figure3(), "figure3_direction.pdf")
