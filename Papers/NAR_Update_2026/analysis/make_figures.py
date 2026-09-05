#!/usr/bin/env python
"""Draft the three main-text figures as one PDF.

Half-page each (7.0 x 4.3 in, double-column width), two panels apiece, covering
the six things the paper has to show: network expansion, structure curation,
pKa integration, reaction-energy generation and comparison, direction-prediction
approaches compared, and atom mapping.

Every number is measured from the database or from this session's audits; none
is illustrative. Sources are named in NUMBERS below so a reviewer can retrace
each one. Where a figure would need a result that does not exist yet -- the
direction-sensitivity study over the model corpus is still \\TBD in the draft --
the panel shows what IS measurable (cross-source direction agreement over the
database) rather than inventing the missing study.

Colours are the validated categorical palette; direction uses a diverging
encoding because the variable is polarity (forward / reversible / reverse), and
the 2020-vs-2026 pairs use two steps of one hue because that is magnitude over
time, not identity.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "figures" / "main_figures_draft.pdf"

# ---- palette (validated: see dataviz references/palette.md) ----------------
BLUE, ORANGE, AQUA, YELLOW, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"
BLUE_200, BLUE_350 = "#9ec5f4", "#5598e7"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8983"
GRID, SURFACE = "#e8e7e3", "#ffffff"
NEUTRAL = "#c9c8c2"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.2,
    "axes.edgecolor": GRID, "axes.linewidth": 0.6,
    "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 6.8, "ytick.labelsize": 6.8,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})

# ---- measured values -------------------------------------------------------
NUMBERS = {
    # MANUSCRIPT.md Table 3 (2020 = Seaver et al., NAR 2021). COMPOUNDS, and
    # only sources that actually supply structures. MetaNetX and BiGG are
    # excluded as mapping / model namespaces rather than structure providers;
    # Rhea is excluded because its 207 InChIKeys alias to zero ModelSEED
    # compounds -- it is a reaction resource here.
    "growth": [("MetaCyc", 19138, 25740), ("KEGG", 17760, 17803),
               ("ChEBI", 0, 11429)],
    # ModelSEED compounds for which each source supplies a structure. These
    # overlap: a compound may carry one from several sources, so they do not
    # sum to the 36,943 total.
    "struct_src": [("MetaCyc", 18801), ("KEGG", 15278), ("ChEBI", 9446)],
    "struct_total": (36943, 45708),
    "coverage": [("2020", 28120, 33992), ("2026", 36943, 45708)],
    # provenance audit of the ADOPTED (collapse-gated) cache, 2026-09-04.
    # ChemAxon-derived pools the cache tier (2,238), our own Marvin table
    # (1,497) and Zenodo carry-over (525): all three descend from cxcalc, and
    # splitting them in the figure implies a distinction the reader cannot use.
    "pka_cpd": [("MolGpKa", 18923, BLUE), ("ChemAxon-derived", 5038, ORANGE),
                ("IUPAC", 1249, AQUA), ("other / unattributed", 19, VIOLET)],
    # distinct scored reactions classified by their compounds' protonation
    # provenance -- NOT compound-reaction incidences, which double-count a
    # reaction once per affected compound
    "pka_traffic": [("all compounds open", 3186, BLUE),
                    ("\u2265 1 ChemAxon-derived compound", 21725, ORANGE),
                    ("no ionizable compound", 159, NEUTRAL)],
    # Biochemistry/*.json thermodynamics dicts, sentinel (1e7) rows EXCLUDED:
    # Group contribution stores an entry for essentially every reaction but
    # 26,555 of them are the 10000000.0 placeholder, so the raw entry count
    # reads as 100% coverage and is not coverage at all.
    "energy_rxn": [("dGPredictor", 29617), ("Group contribution", 29447),
                   ("eQuilibrator", 25070)],
    "energy_total": 56012,
    # direction derived per source from the same dicts
    "direction": [("eQuilibrator", 8650, 10957, 1071),
                  ("Group contribution", 10514, 17510, 1423),
                  ("dGPredictor", 6870, 20723, 2024)],
    "agree": [("all three agree", 10702), ("two of three", 8727), ("all differ", 202)],
    # SourceGrading harness, run 2026-09-04 against the TECRDB anchor
    "dir_grade": [("GOLD", 2185, 2337, 88), ("SILVER", 5156, 8172, 606),
                  ("BRONZE", 5262, 6654, 1230)],
    "dir_src": [("GOLD", 3357, 489, 764), ("SILVER", 3308, 10304, 322),
                ("BRONZE", 160, 10284, 2702)],
    # Structures/AtomMappings/rxns_confidence.tsv
    "atom": [("clean", 25058, BLUE), ("salvaged", 7819, BLUE_350),
             ("not mapped", 56012 - 32877, NEUTRAL)],
}


def strip(ax, keep_x=True):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(keep_x)
    ax.tick_params(length=0)
    ax.set_axisbelow(True)


def panel_tag(ax, letter, title, tx=None):
    """Letter and title. tx must grow on narrow panels or the two collide."""
    ax.text(-0.055, 1.16, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="top", ha="left", color=INK)
    if tx is None:
        w = ax.get_position().width
        tx = 0.022 / max(w, 0.08) * 0.30
    ax.text(tx, 1.16, title, transform=ax.transAxes, fontsize=7.6,
            va="top", ha="left", color=INK)


def swatches(ax, items, y=-0.15, x0=0.0, dx=None, size=6.4, vertical=False):
    """Inline colour key. Identity must never rest on colour alone, and the
    validator's contrast WARN on the lighter hues obliges a visible label."""
    x, yy = x0, y
    for name, col in items:
        w = 0.020 if not vertical else 0.026
        ax.add_patch(plt.Rectangle((x, yy), w, 0.030, transform=ax.transAxes,
                                   facecolor=col, edgecolor="none", clip_on=False,
                                   zorder=5))
        ax.text(x + w + 0.010, yy + 0.014, name, transform=ax.transAxes, fontsize=size,
                va="center", ha="left", color=INK2, clip_on=False)
        if vertical:
            yy -= 0.058
        else:
            x += (dx if dx else 0.034 + 0.0115 * len(name))


def k(n):
    return f"{n/1000:.1f}k" if n >= 1000 else str(n)


# ============================ FIGURE 1 ======================================
def figure1():
    """Structure sources only: who supplies structures, and how far they reach."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 4.3))
    fig.subplots_adjust(left=0.115, right=0.975, top=0.80, bottom=0.14, wspace=0.42)

    rows = NUMBERS["growth"]; ys = range(len(rows))[::-1]; h = 0.32
    for i, (lab, old, new_) in zip(ys, rows):
        if old:
            a.barh(i + h / 2 + 0.02, old, height=h, color=BLUE_200, zorder=3)
            a.text(old + 400, i + h / 2 + 0.02, k(old), va="center", fontsize=6.4, color=MUTED)
            pct = f"+{100*(new_-old)/old:.0f}%"
        else:
            a.text(400, i + h / 2 + 0.02, "not a source in 2020", va="center",
                   fontsize=6.2, color=MUTED, style="italic")
            pct = "new"
        a.barh(i - h / 2 - 0.02, new_, height=h, color=BLUE, zorder=3)
        a.text(new_ + 400, i - h / 2 - 0.02, k(new_), va="center", fontsize=6.4, color=INK)
        a.text(-800, i, lab, va="center", ha="right", fontsize=7.2, color=INK)
        a.text(35500, i, pct, va="center", ha="right", fontsize=6.8,
               color=INK2, fontweight="bold")
    a.set_yticks([]); a.set_xlim(0, 36000); a.set_ylim(-0.7, len(rows) - 0.3)
    a.set_xticks([0, 10000, 20000, 30000]); a.set_xticklabels(["0", "10k", "20k", "30k"])
    a.xaxis.grid(True, color=GRID, lw=0.6); strip(a)
    swatches(a, [("2020", BLUE_200), ("2026", BLUE)], y=1.02, x0=0.60, dx=0.135)
    panel_tag(a, "A", "Compounds per structure source")

    have, tot = NUMBERS["struct_total"]
    rows = NUMBERS["struct_src"]; ys = range(len(rows))[::-1]
    for i, (lab, v) in zip(ys, rows):
        b.barh(i, v, height=0.42, color=BLUE, zorder=3)
        b.text(v + 400, i, f"{k(v)}  {100*v/have:.0f}%", va="center", fontsize=6.6, color=INK)
        b.text(-700, i, lab, va="center", ha="right", fontsize=7.2, color=INK)
    b.axvline(have, color=ORANGE, lw=1.3, zorder=4)
    b.text(have - 500, len(rows) - 0.62, f"{k(have)} compounds\nwith a structure",
           ha="right", va="top", fontsize=6.4, color=ORANGE, fontweight="bold")
    b.set_yticks([]); b.set_xlim(0, 41000); b.set_ylim(-0.7, len(rows) - 0.25)
    b.set_xticks([0, 10000, 20000, 30000]); b.set_xticklabels(["0", "10k", "20k", "30k"])
    b.xaxis.grid(True, color=GRID, lw=0.6); strip(b)
    b.text(0.0, -0.115, f"sources overlap; {100*have/tot:.0f}% of all {k(tot)} compounds "
           f"carry a structure", transform=b.transAxes, fontsize=6.3, color=MUTED)
    panel_tag(b, "B", "Structures reaching ModelSEED compounds")
    return fig


# ============================ FIGURE 2 ======================================
def figure2():
    import json
    fig = plt.figure(figsize=(7.0, 4.3))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 0.75, 1.10],
                          left=0.20, right=0.975, top=0.90, bottom=0.115,
                          hspace=1.05, wspace=0.34)
    a = fig.add_subplot(gs[0, :]); b = fig.add_subplot(gs[1, :])

    for row, (key, label) in enumerate([("pka_traffic", "by scored reaction"),
                                        ("pka_cpd", "by compound")]):
        segs = NUMBERS[key]; tot = sum(v for _, v, _ in segs); x = 0
        for name, v, col in segs:
            pct = v / tot * 100
            a.barh(row, pct, left=x, height=0.42, color=col, zorder=3,
                   edgecolor=SURFACE, lw=1.2)
            if pct > 24:
                a.text(x + pct / 2, row, f"{name}  {pct:.0f}%", ha="center", va="center",
                       fontsize=6.3, color="white", fontweight="bold")
            elif pct > 11:
                a.text(x + pct / 2, row, f"{pct:.0f}%", ha="center", va="center",
                       fontsize=6.3, color="white", fontweight="bold")
            x += pct
        a.text(-1.5, row, label, ha="right", va="center", fontsize=7.0, color=INK)
    a.set_xlim(0, 100); a.set_ylim(-0.55, 1.55); a.set_yticks([])
    a.set_xticks([0, 25, 50, 75, 100]); a.set_xticklabels(["0", "25", "50", "75", "100%"])
    a.xaxis.grid(True, color=GRID, lw=0.6); strip(a)
    swatches(a, [(n, c) for n, _v, c in NUMBERS["pka_cpd"]], y=-0.42, x0=0.0, size=6.1)
    panel_tag(a, "A", "Protonation source")

    rows = NUMBERS["energy_rxn"]; tot = NUMBERS["energy_total"]
    ys = range(len(rows))[::-1]
    for i, (lab, v) in zip(ys, rows):
        b.barh(i, v, height=0.46, color=BLUE, zorder=3)
        b.barh(i, tot - v, left=v + 260, height=0.46, color=NEUTRAL, zorder=3, alpha=0.55)
        b.text(v - 700, i, f"{k(v)}  ({100*v/tot:.0f}%)", va="center", ha="right",
               fontsize=6.6, color="white", fontweight="bold")
        b.text(-900, i, lab, va="center", ha="right", fontsize=7.0, color=INK)
    b.set_yticks([]); b.set_xlim(0, tot * 1.02); b.set_ylim(-0.6, len(rows) - 0.4)
    b.set_xticks([0, 20000, 40000, 56012]); b.set_xticklabels(["0", "20k", "40k", "56k"])
    b.xaxis.grid(True, color=GRID, lw=0.6); strip(b)
    panel_tag(b, "B", "Reactions with a real energy, of 56,012 (sentinels excluded)")

    # C -- reported uncertainty, one axis per source. Deliberately NOT shared:
    # the three report on different scales, and forcing a common axis would
    # flatten two of them into a spike against eQuilibrator's tail.
    sig = json.loads((Path(__file__).resolve().parent
                      / "figure_data_sigma.json").read_text())
    order = [("eQuilibrator", BLUE), ("Group contribution", AQUA), ("dGPredictor", VIOLET)]
    for j, (name, col) in enumerate(order):
        ax = fig.add_subplot(gs[2, j])
        v = sig[name]["vals"]
        hi = sorted(v)[int(0.97 * len(v))]      # clip the tail, then bin INSIDE
        ax.hist([x for x in v if x <= hi], bins=30, range=(0, hi),
                color=col, zorder=3, linewidth=0)
        med = v[len(v) // 2]
        ax.axvline(med, color=INK, lw=0.9, zorder=4)
        ax.text(0.97, 0.90, f"median {med:.1f}", transform=ax.transAxes, ha="right",
                fontsize=6.2, color=INK, fontweight="bold")
        ax.set_title(name, fontsize=6.8, color=INK, pad=3)
        ax.set_yticks([]); ax.tick_params(labelsize=6.0)
        ax.set_xlim(0, hi)
        strip(ax, keep_x=True)
        if sig[name]["dropped"]:
            ax.text(0.0, -0.30, f"{sig[name]['dropped']:,} with no estimate",
                    transform=ax.transAxes, fontsize=5.8, color=MUTED)
        if j == 0:
            ax.text(-0.10, 1.42, "C", transform=ax.transAxes, fontsize=9.5,
                    fontweight="bold", va="top", color=INK)
            ax.text(0.16, 1.42, "Reported uncertainty (kcal/mol), separate axes",
                    transform=ax.transAxes, fontsize=7.6, va="top", color=INK)
    return fig


# ============================ FIGURE 3 ======================================
def figure3():
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 4.3),
                                  gridspec_kw={"width_ratios": [1.15, 1.15, 0.62]})
    fig.subplots_adjust(left=0.095, right=0.975, top=0.80, bottom=0.155, wspace=0.52)

    FWD, REV, BACK = "#1c5cab", NEUTRAL, "#c2410c"
    rows = NUMBERS["dir_grade"]; ys = range(len(rows))[::-1]
    for i, (lab, f, e, r) in zip(ys, rows):
        tot = f + e + r; x = 0
        for v, col, nm in ((f, FWD, "\u2192"), (e, REV, "\u2194"), (r, BACK, "\u2190")):
            pct = v / tot * 100
            a.barh(i, pct, left=x, height=0.44, color=col, zorder=3,
                   edgecolor=SURFACE, lw=1.2)
            if pct > 11:
                a.text(x + pct / 2, i, f"{nm} {pct:.0f}%", ha="center", va="center",
                       fontsize=6.6, fontweight="bold",
                       color="white" if col != NEUTRAL else INK)
            x += pct
        a.text(-1.8, i, lab, va="center", ha="right", fontsize=7.2, color=INK)
        a.text(104.0, i, k(tot), va="center", fontsize=6.4, color=MUTED)
    a.set_xlim(0, 100); a.set_ylim(-0.75, len(rows) - 0.3); a.set_yticks([])
    a.set_xticks([0, 50, 100]); a.set_xticklabels(["0", "50", "100%"])
    a.xaxis.grid(True, color=GRID, lw=0.6); strip(a)
    swatches(a, [("forward", FWD), ("reversible", REV), ("reverse", BACK)],
             y=-0.115, x0=0.0, vertical=True, size=6.2)
    panel_tag(a, "A", "Direction by evidence grade")

    SC = [("eQuilibrator", BLUE), ("Group contribution", AQUA),
          ("dGPredictor", VIOLET)]
    rows = NUMBERS["dir_src"]; ys = range(len(rows))[::-1]
    for i, (lab, *vals) in zip(ys, rows):
        tot = sum(vals); x = 0
        for v, (nm, col) in zip(vals, SC):
            pct = v / tot * 100
            b.barh(i, pct, left=x, height=0.44, color=col, zorder=3,
                   edgecolor=SURFACE, lw=1.2)
            if pct > 13:
                b.text(x + pct / 2, i, f"{pct:.0f}%", ha="center", va="center",
                       fontsize=6.6, color="white", fontweight="bold")
            x += pct
        b.text(-1.8, i, lab, va="center", ha="right", fontsize=7.2, color=INK)
    b.set_xlim(0, 100); b.set_ylim(-0.75, len(rows) - 0.3); b.set_yticks([])
    b.set_xticks([0, 50, 100]); b.set_xticklabels(["0", "50", "100%"])
    b.xaxis.grid(True, color=GRID, lw=0.6); strip(b)
    swatches(b, [(n, c) for n, c in SC], y=-0.115, x0=0.0, vertical=True, size=6.2)
    panel_tag(b, "B", "Which source earned the grade")

    segs = NUMBERS["atom"]; tot = sum(v for _, v, _ in segs); base = 0
    for name, v, col in segs:
        c.bar(0, v, bottom=base + (260 if base else 0), width=0.55, color=col, zorder=3)
        c.text(0.36, base + v / 2, f"{name}\n{k(v)}  {100*v/tot:.0f}%", va="center",
               ha="left", fontsize=6.3, color=INK if col == NEUTRAL else col,
               fontweight="bold")
        base += v + 260
    c.set_xlim(-0.45, 1.75); c.set_ylim(0, tot * 1.03); c.set_xticks([])
    c.set_yticks([0, 20000, 40000, 56012]); c.set_yticklabels(["0", "20k", "40k", "56k"])
    c.yaxis.grid(True, color=GRID, lw=0.6); strip(c, keep_x=False)
    panel_tag(c, "C", "Atom mapping", tx=0.26)
    return fig


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT) as pdf:
        for fn in (figure1, figure2, figure3):
            fig = fn(); pdf.savefig(fig); plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
