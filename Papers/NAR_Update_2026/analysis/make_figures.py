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

import grace_style as grace

OUT = Path(__file__).resolve().parent.parent / "figures" / "main_figures_draft.pdf"

# ---- palette (validated: see dataviz references/palette.md) ----------------
BLUE, ORANGE, AQUA, YELLOW, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"
BLUE_200, BLUE_350 = "#9ec5f4", "#5598e7"
INK, INK2, MUTED = "#000000", "#1a1a1a", "#4d4d4d"
GRID, SURFACE = "#e8e7e3", "#ffffff"
FRAME = grace.FRAME
NEUTRAL = "#c9c8c2"

plt.rcParams.update(grace.RC)

# ---- provenance, computed from released files ------------------------------
def _pka_provenance():
    """Figure 2A's two bars, derived rather than transcribed.

    Reads Biochemistry/Thermodynamics/ProtonationEvidence/pka_provenance.tsv
    (written by Scripts/Thermodynamics/ProtonationEvidence/build_pka_provenance.py)
    and ModelSEED_Reaction_Energies.tsv. Both ship, so a reader can reproduce
    every percentage in the protonation paragraph.

    REWRITTEN 2026-09-08 for the Marvin 26.1 rebuild. The split is no longer
    open-vs-proprietary -- every resolved ladder is ChemAxon-derived now -- but
    NEW-RUN vs CARRIED-OVER, which is the distinction that survived the change.
    The rebuild reached 67.8% of compounds and 6.1% of scored reactions, because
    the rows it could not rebuild are the cofactors that appear everywhere.
    """
    import csv as _csv, collections as _c, re as _re
    root = Path(__file__).resolve().parents[3]
    prov = root / "Biochemistry/Thermodynamics/ProtonationEvidence/pka_provenance.tsv"
    rxns = root / "Biochemistry/Thermodynamics/eQuilibrator/ModelSEED_Reaction_Energies.tsv"
    rows = [r for r in _csv.DictReader(
        (l for l in prov.open() if not l.startswith("#")), delimiter="\t")]
    eff = _c.Counter(r["effective_source"] for r in rows)
    by_cpd = [("Marvin 26.1 (this release)", eff["marvin"], BLUE),
              ("carried over", eff["carried_over"], ORANGE),
              ("unresolved", eff["unresolved"] + eff["none"], NEUTRAL)]
    src = {r["seed_id"]: r["effective_source"] for r in rows}
    CPD = _re.compile(r"cpd\d{5}")
    a = b = c = 0
    for r in _csv.DictReader((l for l in rxns.open() if not l.startswith("#")), delimiter="\t"):
        if r["status"] != "ok":
            continue
        s = {src.get(x) for x in set(CPD.findall(r["formula"] or ""))}
        if "carried_over" in s: a += 1
        elif "marvin" in s: b += 1
        else: c += 1
    by_rxn = [("all Marvin 26.1", b, BLUE),
              ("\u2265 1 carried over", a, ORANGE),
              ("no resolved ladder", c, NEUTRAL)]
    return by_cpd, by_rxn


def _silver_sigma_band():
    """Figure 2C shading: the sigma span of reactions graded SILVER, per source.

    Reaction-level grade (best_grade in source_grades_wide.tsv) joined back to
    each source's own reported sigma. Answers "what uncertainty does a silver
    reaction actually carry?" -- and shows that the answer is only meaningful
    for two of the three sources. eQuilibrator and dGPredictor separate their
    tiers by sigma (medians 0.26/0.62/1.43 and 1.54/16.49/21.85 for
    gold/silver/bronze); Group contribution does not (9.19/8.99/10.35), which is
    the same flat-error-curve problem that stopped raw sigma being used to rank
    sources in the first place. The band is p5-p95.
    """
    import csv as _csv, collections as _c
    root = Path(__file__).resolve().parents[3]
    g = root / "Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades"
    best = {}
    with (g / "source_grades_wide.tsv").open() as fh:
        for x in _csv.DictReader(fh, delimiter="\t"):
            if x["best_grade"]:
                best[x["rxn"]] = x["best_grade"]
    sig = _c.defaultdict(list)
    with (g / "source_grades.tsv").open() as fh:
        for x in _csv.DictReader(fh, delimiter="\t"):
            if x["source"] == "TECRDB" or best.get(x["rxn"]) != "SILVER":
                continue
            try:
                sig[x["source"]].append(abs(float(x["sigma"])))
            except (TypeError, ValueError):
                pass
    out = {}
    for k, v in sig.items():
        v.sort()
        if len(v) >= 20:
            out[k] = (v[int(0.05 * len(v))], v[int(0.95 * len(v))])
    return out


def _grade_breakdown():
    """Figure 3C: what underpins each evidence grade, on two axes.

    Derived from the shipped Biochemistry/reaction_*.json `thermo-evidence`
    block, so it stays in step with the release rather than the grading run.
    Returns (by_assessment, by_cross) -- each {grade: [(label, n, colour)]}.
    """
    import json as _json, glob as _glob, collections as _c
    root = Path(__file__).resolve().parents[3]
    A = _c.defaultdict(_c.Counter); X = _c.defaultdict(_c.Counter)
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _json.load(open(f)):
            e = r.get("thermo-evidence")
            if not e:
                continue
            A[e["grade"]][e["assessment"]] += 1
            X[e["grade"]][e.get("cross-source", "no cross-check")] += 1
    AC = {"measured": ORANGE, "self-certain": BLUE, "self-confident": BLUE_350,
          "unconfident": NEUTRAL}
    XC = {"corroborated": BLUE, "outvoted": ORANGE, "unpaired": AQUA,
          "no cross-check": NEUTRAL}
    order_a = ["measured", "self-certain", "self-confident", "unconfident"]
    order_x = ["corroborated", "outvoted", "unpaired", "no cross-check"]
    by_a = {g: [(k, A[g][k], AC[k]) for k in order_a if A[g][k]] for g in A}
    by_x = {g: [(k, X[g][k], XC[k]) for k in order_x if X[g][k]] for g in X}
    return by_a, by_x


def _direction_counts():
    """Figure 3A: direction assigned by each source, DERIVED from the shipped
    reaction JSON. Four states, not three -- "undetermined" is the one the 2020
    release folded into "reversible" and is the point of the panel."""
    import json as _json, glob as _glob, collections as _c
    root = Path(__file__).resolve().parents[3]
    counts = {s: _c.Counter() for s in
              ("eQuilibrator", "Group contribution", "dGPredictor")}
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _json.load(open(f)):
            for s, v in (r.get("thermodynamics") or {}).items():
                if s in counts and isinstance(v, list) and len(v) > 2:
                    counts[s][v[2]] += 1
    return [(s, counts[s][">"], counts[s]["="], counts[s]["<"], counts[s]["?"])
            for s in ("eQuilibrator", "dGPredictor", "Group contribution")]


def _agreement_counts():
    """Figure 3B: eQuilibrator against dGPredictor on the reactions both score."""
    import json as _json, glob as _glob
    root = Path(__file__).resolve().parents[3]
    agree = opp = eq_only = dg_only = neither = partial = 0
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _json.load(open(f)):
            th = r.get("thermodynamics") or {}
            x, y = th.get("eQuilibrator"), th.get("dGPredictor")
            if not (isinstance(x, list) and len(x) > 2
                    and isinstance(y, list) and len(y) > 2):
                continue
            e, d = x[2], y[2]
            if e == "?" and d == "?": neither += 1
            elif d == "?": eq_only += 1
            elif e == "?": dg_only += 1
            elif e == d: agree += 1
            elif {e, d} == {">", "<"}: opp += 1
            else: partial += 1     # one calls a direction, the other reversible
    # Six mutually exclusive states that must sum to the shared total. The
    # "partial" bucket is 654 reactions and was silently dropped in the first
    # draft of this function, which made the panel sum to 97.4%.
    return [("both agree", agree, BLUE), ("only eQuilibrator", eq_only, BLUE_350),
            ("only dGPredictor", dg_only, AQUA), ("neither", neither, NEUTRAL),
            ("one direction, one reversible", partial, YELLOW),
            ("opposite", opp, ORANGE)]


# ---- measured values -------------------------------------------------------
NUMBERS = {
    # MANUSCRIPT.md Table 3 (untracked; local to the author's tree, not in the repository)
    # (2020 = Seaver et al., NAR 2021). COMPOUNDS, and
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
    # Protonation provenance. READ FROM THE SHIPPED TABLE, not transcribed --
    # see _pka_provenance() below. These were hardcoded constants until
    # 2026-09-07 and were WRONG: they counted the pKa cascade's answers and
    # treated carry-over as 525 compounds when the cache holds 4,122 rows from
    # the pinned release. Correcting it moves the layer from a reported 20%
    # ChemAxon-derived to 28%, and 87% of scored reactions to 95% -- the
    # proprietary dependency is LARGER than the paper claimed. Task #58.
    "pka_cpd": None,      # filled by _pka_provenance()
    "pka_traffic": None,  # filled by _pka_provenance()
    # Biochemistry/*.json thermodynamics dicts, sentinel (1e7) rows EXCLUDED:
    # Group contribution stores an entry for essentially every reaction but
    # 26,555 of them are the 10000000.0 placeholder, so the raw entry count
    # reads as 100% coverage and is not coverage at all.
    "energy_rxn": [("dGPredictor", 29617), ("Group contribution", 29447),
                   ("eQuilibrator", 21789)],
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


NUMBERS["pka_cpd"], NUMBERS["pka_traffic"] = _pka_provenance()
NUMBERS["direction"] = _direction_counts()
NUMBERS["agreement"] = _agreement_counts()
NUMBERS["silver_sigma"] = _silver_sigma_band()
NUMBERS["grade_assess"], NUMBERS["grade_cross"] = _grade_breakdown()


def strip(ax, keep_x=True, value_axis="x"):
    """Grace frame. Named strip() still because every call site says strip();
    it no longer strips anything -- Grace closes the box rather than opening it.
    keep_x is vestigial and ignored: a Grace axes always has all four sides."""
    return grace.frame(ax, value_axis=value_axis)


def panel_tag(ax, letter, title, tx=None):
    """Letter and title, offset in POINTS above the frame.

    Previously offset by axes fraction, which is height-dependent: the same
    0.16 that cleared a tall panel put figure 2B's title on top of its frame.
    tx is retained for call-site compatibility and ignored.
    """
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(-7, 8), textcoords="offset points",
                fontsize=9.5, fontweight="bold", va="baseline", ha="left",
                color=INK, annotation_clip=False)
    ax.annotate(title, xy=(0, 1), xycoords="axes fraction",
                xytext=(8, 8), textcoords="offset points",
                fontsize=7.6, va="baseline", ha="left", color=INK,
                annotation_clip=False)


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

    rows = NUMBERS["growth"]; ys = range(len(rows))[::-1]; h = 0.38
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
    strip(a)
    # Was above the frame at y=1.02, where it collided with the panel title once
    # panel_tag moved to a fixed point offset. Below the axis, as in figure 3.
    swatches(a, [("2020", BLUE_200), ("2026", BLUE)], y=-0.105, x0=0.0, dx=0.135)
    panel_tag(a, "A", "Compounds per structure source")

    have, tot = NUMBERS["struct_total"]
    rows = NUMBERS["struct_src"]; ys = range(len(rows))[::-1]
    for i, (lab, v) in zip(ys, rows):
        b.barh(i, v, height=0.56, color=BLUE, zorder=3)
        b.text(v + 400, i, f"{k(v)}  {100*v/have:.0f}%", va="center", fontsize=6.6, color=INK)
        b.text(-700, i, lab, va="center", ha="right", fontsize=7.2, color=INK)
    b.axvline(have, color=ORANGE, lw=1.3, zorder=4)
    b.text(have - 500, len(rows) - 0.62, f"{k(have)} compounds\nwith a structure",
           ha="right", va="top", fontsize=6.4, color=ORANGE, fontweight="bold")
    b.set_yticks([]); b.set_xlim(0, 41000); b.set_ylim(-0.7, len(rows) - 0.25)
    b.set_xticks([0, 10000, 20000, 30000]); b.set_xticklabels(["0", "10k", "20k", "30k"])
    strip(b)
    # Two lines, not one: as a single line this ran past the figure edge and
    # was clipped mid-word in the typeset PDF.
    b.text(0.0, -0.105, f"sources overlap; {100*have/tot:.0f}% of all\n"
           f"{k(tot)} compounds carry a structure",
           transform=b.transAxes, fontsize=6.3, color=MUTED,
           va="top", linespacing=1.35)
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
            a.barh(row, pct, left=x, height=0.54, color=col, zorder=3,
                   edgecolor=FRAME, lw=0.45)
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
    strip(a)
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
    strip(b)
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
        # Shade the sigma span of SILVER-graded reactions, so a reader can see
        # what uncertainty a tier actually corresponds to -- and, for group
        # contribution, that it corresponds to almost the whole distribution.
        band = NUMBERS["silver_sigma"].get(name)
        if band:
            lo, bhi = band
            ax.axvspan(lo, min(bhi, hi), color=YELLOW, alpha=0.20, lw=0, zorder=2)
            for xv in (lo, bhi):
                if xv <= hi:
                    ax.axvline(xv, color=YELLOW, lw=0.7, zorder=2.5)
        med = v[len(v) // 2]
        ax.axvline(med, color=INK, lw=0.9, zorder=4)
        # White bbox: the median rule runs up through this label's line.
        ax.text(0.97, 0.86, f"median {med:.1f}", transform=ax.transAxes, ha="right",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.8),
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
            # Figure-level so it cannot collide with the per-panel titles.
            fig.text(0.205, 0.016, "shaded band: p5\u2013p95 of the reactions this "
                     "source grades silver", fontsize=6.1, color=MUTED)

    return fig


# ============================ FIGURE 3 ======================================
def figure3():
    fig = plt.figure(figsize=(7.0, 5.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 0.62],
                          height_ratios=[1.0, 0.88],
                          left=0.135, right=0.945, top=0.90, bottom=0.085,
                          wspace=0.80, hspace=0.62)
    a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[0, 2]); d = fig.add_subplot(gs[1, :])

    FWD, REV, BACK = "#1c5cab", NEUTRAL, "#c2410c"
    # A -- direction assigned by each SOURCE independently (was: by grade).
    # FOUR states, not three. "undetermined" is the one the 2020 release folded
    # into "reversible", and separating them is the point of the panel: 63% of
    # dGPredictor's calls live there.
    UND = NEUTRAL
    rows = NUMBERS["direction"]; ys = range(len(rows))[::-1]
    for i, (lab, f, e, r, q) in zip(ys, rows):
        tot = f + e + r + q; x = 0
        for v, col, nm in ((f, FWD, "\u2192"), (e, AQUA, "\u2194"),
                           (r, BACK, "\u2190"), (q, UND, "?")):
            pct = v / tot * 100
            a.barh(i, pct, left=x, height=0.60, color=col, zorder=3,
                   edgecolor=FRAME, lw=0.45)
            if pct > 9:
                a.text(x + pct / 2, i, f"{nm} {pct:.0f}%", ha="center", va="center",
                       fontsize=6.4, fontweight="bold",
                       color="white" if col != UND else INK)
            x += pct
        a.text(-3.5, i, lab.replace("Group contribution", "Group contrib."),
               va="center", ha="right", fontsize=7.2, color=INK)
        a.text(102.5, i, k(tot), va="center", ha="left", fontsize=6.2, color=MUTED)
    a.set_xlim(0, 100); a.set_ylim(-0.62, len(rows) - 0.38); a.set_yticks([])
    a.set_xticks([0, 50, 100]); a.set_xticklabels(["0", "50", "100%"])
    strip(a)
    swatches(a, [("forward", FWD), ("reversible", AQUA), ("reverse", BACK),
                 ("undetermined", UND)], y=-0.115, x0=0.0, vertical=True, size=6.2)
    panel_tag(a, "A", "Direction by source")

    # B -- eQuilibrator against dGPredictor on the 24,804 reactions both score.
    # Ranked bars rather than one stacked bar: the categories span 13,575 to 198
    # and the small ones are the interesting ones. Six states, summing to the
    # shared total exactly.
    SHORT = {"one direction, one reversible": "one dir., one rev."}
    segs = sorted(NUMBERS["agreement"], key=lambda s: -s[1])
    tot = sum(v for _, v, _ in segs)
    ys = range(len(segs))[::-1]
    for i, (nm, v, col) in zip(ys, segs):
        b.barh(i, v, height=0.62, color=col, zorder=3, edgecolor=FRAME, lw=0.45)
        b.text(v + tot * 0.015, i, f"{k(v)}  {100*v/tot:.1f}%", va="center",
               ha="left", fontsize=6.2, color=INK)
        b.text(-tot * 0.02, i, SHORT.get(nm, nm), va="center", ha="right",
               fontsize=6.6, color=INK)
    b.set_xlim(0, tot * 0.78); b.set_ylim(-0.62, len(segs) - 0.38); b.set_yticks([])
    b.set_xticks([0, 5000, 10000, 15000]); b.set_xticklabels(["0", "5k", "10k", "15k"])
    strip(b)
    panel_tag(b, "B", "eQuilibrator vs dGPredictor")

    segs = NUMBERS["atom"]; tot = sum(v for _, v, _ in segs); base = 0
    for name, v, col in segs:
        c.bar(0, v, bottom=base + (260 if base else 0), width=0.55, color=col, zorder=3)
        c.text(0.36, base + v / 2, f"{name}\n{k(v)}  {100*v/tot:.0f}%", va="center",
               ha="left", fontsize=6.3, color=INK if col == NEUTRAL else col,
               fontweight="bold")
        base += v + 260
    c.set_xlim(-0.45, 1.75); c.set_ylim(0, tot * 1.03); c.set_xticks([])
    c.set_yticks([0, 20000, 40000, 56012]); c.set_yticklabels(["0", "20k", "40k", "56k"])
    strip(c, keep_x=False, value_axis="y")
    panel_tag(c, "C", "Atom mapping")

    # D -- what underpins each grade, on two axes. Six lanes: each grade split
    # by the deciding source's self-assessment, then the same three split by
    # what the other sources made of it. Read together they say where a tier's
    # authority comes from -- gold from confidence plus corroboration, bronze
    # from neither.
    GR = ["gold", "silver", "bronze"]
    lanes = [(g, NUMBERS["grade_assess"].get(g, [])) for g in GR] + \
            [(g, NUMBERS["grade_cross"].get(g, [])) for g in GR]
    ys = range(len(lanes))[::-1]
    for i, (lab, segs) in zip(ys, lanes):
        tot = sum(v for _, v, _ in segs) or 1
        x = 0
        for nm, v, col in segs:
            pct = v / tot * 100
            d.barh(i, pct, left=x, height=0.60, color=col, zorder=3,
                   edgecolor=FRAME, lw=0.45)
            if pct > 14:
                d.text(x + pct / 2, i, f"{nm} {pct:.0f}%", ha="center", va="center",
                       fontsize=5.9, fontweight="bold",
                       color="white" if col != NEUTRAL else INK)
            x += pct
        d.text(-1.8, i, lab, va="center", ha="right", fontsize=7.0, color=INK)
        d.text(102.5, i, k(tot), va="center", ha="left", fontsize=6.0, color=MUTED)
    d.set_xlim(0, 100); d.set_ylim(-0.62, len(lanes) - 0.38); d.set_yticks([])
    d.set_xticks([0, 50, 100]); d.set_xticklabels(["0", "50", "100%"])
    strip(d)
    d.text(-0.118, 0.80, "by self-assessment", transform=d.transAxes,
           rotation=90, va="center", ha="center", fontsize=6.1, color=MUTED)
    d.text(-0.118, 0.26, "by cross-source", transform=d.transAxes,
           rotation=90, va="center", ha="center", fontsize=6.1, color=MUTED)
    panel_tag(d, "D", "What underpins each evidence grade")
    return fig


def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT) as pdf:
        for fn in (figure1, figure2, figure3):
            fig = fn(); pdf.savefig(fig); plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
