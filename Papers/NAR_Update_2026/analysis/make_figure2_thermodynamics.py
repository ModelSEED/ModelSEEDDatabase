#!/usr/bin/env python
"""Figure 2: provenance, coverage and reported uncertainty of the thermodynamic sources.

Regenerate with:  python make_figure2_thermodynamics.py
Shared palette, NUMBERS and helpers live in figure_common.py.
"""
from figure_common import *  # noqa: F401,F403 -- palette, NUMBERS, helpers
from figure_common import save


# ============================ FIGURE 2 ======================================
def figure2():
    """Two columns: the two coverage bars on the left, the three uncertainty
    distributions stacked as a column on the right. Previously three full-width
    rows, which made the figure tall; side by side it spans the page instead and
    costs roughly a quarter of the vertical space."""
    import json
    fig = plt.figure(figsize=(7.0, 3.05))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.78, 1.0],
                             left=0.093, right=0.988, top=0.986, bottom=0.088,
                             wspace=0.105)
    left = outer[0, 0].subgridspec(2, 1, hspace=0.318, height_ratios=[1.0, 1.34])
    right = outer[0, 1].subgridspec(3, 1, hspace=0.392)
    a = fig.add_subplot(left[0]); b = fig.add_subplot(left[1])

    def tag(ax, letter):
        ax.annotate(letter, xy=(0.98, 0.94), xycoords="axes fraction", fontsize=9.0,
                    fontweight="bold", va="top", ha="right", color=INK, zorder=6)

    # counts axis: both rows on one scale, so reactions and compounds are
    # directly comparable -- the thing percentages hid (Sam 2026-09-14)
    AMAX = max(sum(s[1] for s in NUMBERS[k_]) for k_ in ("ladder_rxn", "ladder_cpd"))
    for row, (key, label) in enumerate([("ladder_rxn", "reactions"),
                                        ("ladder_cpd", "compounds")]):
        segs = NUMBERS[key]; tot = sum(s[1] for s in segs); x = 0
        for name, v, col, hatch in segs:
            pct = v                      # counts, not percent
            # hatch marks the SMILES-derived route through the same Marvin
            # release; white strokes read against both the blue and the orange
            a.barh(row, pct, left=x, height=0.74, color=col, zorder=3,
                   edgecolor=SURFACE if hatch else FRAME,
                   lw=0.45, hatch=hatch or None)
            # ink on the pale fills, white on the saturated ones; a tight
            # surface-coloured box lifts the number clear of the hatch strokes
            fg = INK if col in (GRID, NEUTRAL) else "white"
            box = dict(facecolor=col, edgecolor="none", pad=0.9) if hatch else None
            if pct > AMAX * 0.22:
                a.text(x + pct / 2, row, f"{name}  {k(v)}", ha="center", va="center",
                       fontsize=6.0, color=fg, fontweight="bold", zorder=5, bbox=box)
            elif pct > AMAX * 0.09:
                a.text(x + pct / 2, row, k(v), ha="center", va="center",
                       fontsize=6.0, color=fg, fontweight="bold", zorder=5, bbox=box)
            x += pct
        a.text(-AMAX * 0.015, row, label, ha="right", va="center", fontsize=7.0, color=INK)
    a.set_xlim(0, AMAX); a.set_ylim(-0.62, 2.42); a.set_yticks([])
    _end = NUMBERS["energy_total"]
    a.set_xticks([0, 20000, 40000, _end])
    a.set_xticklabels(["0", "20k", "40k", f"{_end//1000}k"])
    strip(a)
    # Four segments, two of them previously identified only in the caption.
    # The key names all four; hatch is reproduced in its own swatch so the
    # SMILES route is not left to be inferred from the fill.
    from matplotlib.patches import Patch
    a.legend(handles=[
        Patch(facecolor=BLUE, edgecolor=FRAME, lw=0.45, label="Marvin"),
        Patch(facecolor=BLUE, edgecolor=SURFACE, lw=0.45, hatch="xxx",
              label="Marvin, from SMILES"),
        Patch(facecolor=ORANGE, edgecolor=FRAME, lw=0.45, label="no ionizable site"),
        Patch(facecolor=GRID, edgecolor=FRAME, lw=0.45, label="no structure")],
        loc="upper left", frameon=False, fontsize=5.9, ncol=4, handlelength=1.1,
        handleheight=0.85, handletextpad=0.4, columnspacing=0.9, borderpad=0.1,
        borderaxespad=0.15, labelcolor=INK2)
    tag(a, "A")

    rows = NUMBERS["energy_rxn"]; tot = NUMBERS["energy_total"]
    ys = range(len(rows))[::-1]
    SHORT = {"Group contribution": "Group contr."}
    for i, (lab, v) in zip(ys, rows):
        b.barh(i, v, height=0.70, color=BLUE, zorder=3)
        # grey remainder removed 2026-09-14: the count axis already shows the
        # shortfall against the total, so the bar drew the same fact twice
        b.text(v - 700, i, f"{100*v/tot:.0f}%", va="center", ha="right",
               fontsize=6.4, color="white", fontweight="bold")
        b.text(-900, i, SHORT.get(lab, lab), va="center", ha="right",
               fontsize=7.0, color=INK)
    b.set_yticks([]); b.set_xlim(0, tot * 1.02); b.set_ylim(-0.62, len(rows) - 0.38)
    b.set_xticks([0, 20000, 40000, _end])
    b.set_xticklabels(["0", "20k", "40k", f"{_end//1000}k"])
    strip(b)
    tag(b, "B")

    # C -- reported uncertainty, one axis per source, stacked as a right-hand
    # column. Axes are deliberately NOT shared: the three report on different
    # scales, and a common axis would flatten two into a spike against
    # eQuilibrator's tail. Source names sit inside the frame; three stacked
    # titles would cost more height than the panels themselves.
    sig = json.loads((Path(__file__).resolve().parent
                      / "figure_data_sigma.json").read_text())
    order = [("eQuilibrator", BLUE), ("Group contribution", AQUA), ("dGPredictor", VIOLET)]
    SHORT_C = {}
    for j, (name, col) in enumerate(order):
        ax = fig.add_subplot(right[j])
        v = sig[name]["vals"]
        hi = sorted(v)[int(0.97 * len(v))]      # clip the tail, then bin INSIDE
        ax.hist([x for x in v if x <= hi], bins=30, range=(0, hi),
                color=col, zorder=3, linewidth=0)
        band = NUMBERS["silver_sigma"].get(name)
        if band:
            lo, bhi = band
            ax.axvspan(lo, min(bhi, hi), color=YELLOW, alpha=0.20, lw=0, zorder=2)
            for xv in (lo, bhi):
                if xv <= hi:
                    ax.axvline(xv, color=YELLOW, lw=0.7, zorder=2.5)
        med = v[len(v) // 2]
        ax.axvline(med, color=INK, lw=0.9, zorder=4)
        # all three source names share one left edge; the "C" tag sits at the
        # top RIGHT (0.98), so the first panel never needed to yield this corner
        ax.text(0.030, 0.90, SHORT_C.get(name, name),
                transform=ax.transAxes, ha="left", va="top", fontsize=6.3, color=INK,
                fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", pad=0.6),
                zorder=5)
        ax.set_yticks([]); ax.tick_params(labelsize=5.8, pad=1.5)
        ax.set_xlim(0, hi)
        strip(ax, keep_x=True)
        if j == 0:
            tag(ax, "C")
    return fig


if __name__ == "__main__":
    save(figure2(), "figure2_thermodynamics.pdf")
