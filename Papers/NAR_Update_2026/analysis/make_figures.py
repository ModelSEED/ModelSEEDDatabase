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


def _ladder_vintage():
    """Figure 2A: the full extent of the shipped pKa ladders, and how each was
    reached.

    Marvin 26.1 computed most ladders from the per-source InChI files through
    the cxcalc CLI. That CLI refuses polymers and organometallics as "query
    molecules", so those are computed instead from SMILES through the Java
    PkaPlugin (PR #292) -- the same engine and release, a different route in.
    The SMILES-derived share is hatched rather than given its own colour,
    because it is a provenance distinction inside one release, not a separate
    source: engine parity was measured at median |delta| 0.0000 over 300
    compounds both routes can do.

    Denominator is EVERY compound and EVERY reaction, not the InChI-bearing
    subset the panel used before the SMILES route existed -- the point of the
    panel is now coverage of the database rather than of one structure file.

    Resolution goes through BiochemPy.loadPerSourcePkas and the same
    KEGG > MetaCyc > ChEBI > Rhea cascade Update_Compound_pKas.py applies, so
    the panel cannot drift from what the database ships.
    """
    import csv as _csv, glob as _glob, json as _json, collections as _c, sys as _sys
    root = Path(__file__).resolve().parents[3]
    _sys.path.insert(0, str(root / "Libs" / "Python"))
    from BiochemPy import Compounds
    DBS = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]

    def norm(db, e):
        if db == "ChEBI" and e.startswith("CHEBI_"): return e[len("CHEBI_"):]
        if db == "Rhea" and e.startswith("POLYMER_"): return "POLYMER:" + e[len("POLYMER_"):]
        return e

    ps = Compounds().loadPerSourcePkas(DBS)
    inchi = {db: {norm(db, r["external_id"]) for r in _csv.DictReader(
                 open(root / "Biochemistry/Structures" / db / "inchi.tsv"), dialect="excel-tab")}
             for db in DBS}
    ver = {}
    for db in DBS:
        for f in sorted(_glob.glob(str(root / "Biochemistry/Structures" / db / "pkas/*.tsv"))):
            for r in _csv.DictReader(open(f), dialect="excel-tab"):
                ee, kk = r.get("external_id"), r.get("kind")
                if ee and kk:
                    ver[(db, norm(db, ee), kk, r.get("value"))] = r.get("tool_version")
    al = _c.defaultdict(lambda: _c.defaultdict(list))
    with (root / "Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt").open() as fh:
        rd = _csv.reader(fh, delimiter="\t"); next(rd)
        for r in rd:
            if len(r) >= 3: al[r[0]][r[2]].append(r[1])

    allc, struct = set(), set()
    for f in sorted(_glob.glob(str(root / "Biochemistry/compound_*.json"))):
        for c in _json.load(open(f)):
            allc.add(c["id"])
            if c.get("smiles") or c.get("inchikey"): struct.add(c["id"])

    cls = {}
    for cpd in allc:
        got = None
        for db in DBS:
            for a in al[cpd].get(db, []):
                if (db, a) in ps:
                    vs = {ver.get((db, a, k, v)) for k, v in ps[(db, a)].items()}; vs.discard(None)
                    got = ("23.4" if vs == {"23.4"} else "26.1",
                           "InChI" if a in inchi[db] else "SMILES")
                    break
            if got: break
        cls[cpd] = " ".join(got) if got else ("no site" if cpd in struct else "no structure")

    cc = _c.Counter(cls.values())
    # Grouped by ROUTE, not by release: the 23.4 residue is 400 compounds, all
    # of them SMILES-only ids, so it folds into the hatched share rather than
    # earning a fourth colour. The release split is a sentence in the text.
    by_cpd = [("Marvin", cc["26.1 InChI"] + cc["23.4 InChI"], BLUE, ""),
              ("from SMILES", cc["26.1 SMILES"] + cc["23.4 SMILES"], BLUE, "xxx"),
              ("no ionizable site", cc["no site"], NEUTRAL, ""),
              ("no structure", cc["no structure"], GRID, "")]

    rx = _c.Counter()
    for f in sorted(_glob.glob(str(root / "Biochemistry/reaction_*.json"))):
        for r in _json.load(open(f)):
            st = r.get("stoichiometry")
            if not isinstance(st, list) or not st: continue
            v = {cls.get(x["compound"], "no structure") for x in st}
            if "no structure" in v: rx["inc"] += 1
            elif any(z.endswith("SMILES") for z in v): rx["smi"] += 1
            elif any(z.endswith("InChI") for z in v): rx["inchi"] += 1
            else: rx["none"] += 1
    by_rxn = [("Marvin", rx["inchi"], BLUE, ""),
              ("from SMILES", rx["smi"], BLUE, "xxx"),
              ("no ionizable site", rx["none"], NEUTRAL, ""),
              ("no structure", rx["inc"], GRID, "")]
    return by_cpd, by_rxn


def _growth():
    """Figure 1A: compounds carrying an alias from each structure source,
    2020 against 2026. DERIVED from the shipped alias file on both sides.

    The 2020 column used to be transcribed from an untracked MANUSCRIPT.md and
    was wrong for both sources that existed then -- MetaCyc 19,138 against a
    true 19,172, KEGG 17,760 against 17,793. It is now read from
    Unique_ModelSEED_Compound_Aliases.txt at the last commit of 2020
    (fd6c7849, 2020-11-10), six weeks after the paper appeared online, so the
    figure cannot drift from the repository again.

    Rhea is absent by construction: its 207 InChIKeys alias to zero ModelSEED
    compounds, so it contributes reactions, not structures.
    """
    import csv as _csv, io as _io, subprocess as _sp, collections as _c
    root = Path(__file__).resolve().parents[3]
    C2020 = "fd6c7849891ef4bbeb6eac072f5a6f7adff05b0e"
    REL = "Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt"

    def tally(text):
        seen = _c.defaultdict(set)
        rd = _csv.reader(_io.StringIO(text), delimiter="\t"); next(rd)
        for r in rd:
            if len(r) >= 3:
                seen[r[2]].add(r[0])
        return seen

    old = tally(_sp.run(["git", "show", f"{C2020}:{REL}"], cwd=str(root),
                        capture_output=True, text=True, check=True).stdout)
    new = tally((root / REL).read_text())
    return [(s, len(old[s]), len(new[s])) for s in ("MetaCyc", "KEGG", "ChEBI")]


def _reaction_sources():
    """Figure 1C/D: what each PRIMARY database contributes in reactions.

    A reaction is credited to a source if that source carries an alias for it,
    and is UNIQUE to a source when no other primary database does -- that is the
    set the database would not have without it, which is the question Rhea's
    integration raises. "Complete" means every compound in the stoichiometry
    carries a structure, so the reaction can be mass-balanced and decomposed;
    Rhea identifies its compounds through ChEBI, so its completeness is the
    reach of ChEBI structures into ModelSEED.

    Derived from Biochemistry/Aliases/ and the compound records, both shipped.
    """
    import csv as _csv, collections as _c, glob as _glob, json as _json
    root = Path(__file__).resolve().parents[3]
    struct = set()
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "compound_*.json"))):
        for c in _json.load(open(f)):
            if c.get("smiles") or c.get("inchikey"):
                struct.add(c["id"])
    PRIMARY = ("KEGG", "MetaCyc", "Rhea")
    src = _c.defaultdict(set)
    with (root / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt").open() as fh:
        rd = _csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) >= 3 and row[2] in PRIMARY:
                src[row[0]].add(row[2])
    rx = {}
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for x in _json.load(open(f)):
            rx[x["id"]] = {s["compound"] for s in (x.get("stoichiometry") or [])}
    tot = _c.Counter(); uniq = _c.Counter(); ucomp = _c.Counter()
    for rid, ss in src.items():
        cs = rx.get(rid)
        if cs is None:
            continue
        ok = bool(cs) and cs <= struct
        for s in ss:
            tot[s] += 1
        if len(ss) == 1:
            s = next(iter(ss)); uniq[s] += 1
            if ok: ucomp[s] += 1
    return [(s, tot[s], uniq[s], ucomp[s]) for s in PRIMARY]


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
            X[e["grade"]][e.get("cross-source", "neither way")] += 1
    # Four DISTINCT hues per panel (2026-09-15) -- self-certain/self-confident
    # were two blues and read as one category. NEUTRAL is deliberately shared
    # between the panels: it is the null case in both. BLUE is the only hue
    # reused with different meanings, and the panels carry separate keys.
    AC = {"measured": ORANGE, "self-certain": BLUE, "self-confident": AQUA,
          "unconfident": NEUTRAL}
    XC = {"corroborated": AQUA, "disputed": YELLOW, "unpaired": BLUE,
          "neither way": NEUTRAL}
    # Both orders are Table 1's, so figure and table rank the categories alike.
    order_a = ["measured", "self-certain", "self-confident", "unconfident"]
    order_x = ["corroborated", "disputed", "unpaired", "neither way"]
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
              ("eQuilibrator", "Group contribution", "dGPredictor", "LLMs")}
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _json.load(open(f)):
            for s, v in (r.get("thermodynamics") or {}).items():
                if s in counts and isinstance(v, list) and len(v) > 2:
                    counts[s][v[2]] += 1
    # LLMs last: it carries a direction and no energy, so it sits apart from
    # the three numeric sources and takes no part in the evidence grading.
    return [(s, counts[s][">"], counts[s]["="], counts[s]["<"], counts[s]["?"])
            for s in ("eQuilibrator", "dGPredictor", "Group contribution", "LLMs")]


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
    "growth": None,       # filled by _growth()
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
NUMBERS["growth"] = _growth()
NUMBERS["rxn_sources"] = _reaction_sources()
NUMBERS["ladder_cpd"], NUMBERS["ladder_rxn"] = _ladder_vintage()


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
    a.set_xlim(0, AMAX); a.set_ylim(-0.62, 1.62); a.set_yticks([])
    a.set_xticks([0, 20000, 40000, 56002]); a.set_xticklabels(["0", "20k", "40k", "56k"])
    strip(a)
    tag(a, "A")

    rows = NUMBERS["energy_rxn"]; tot = NUMBERS["energy_total"]
    ys = range(len(rows))[::-1]
    SHORT = {"Group contribution": "Group contr."}
    for i, (lab, v) in zip(ys, rows):
        b.barh(i, v, height=0.70, color=BLUE, zorder=3)
        # grey remainder removed 2026-09-14: the count axis already shows the
        # shortfall against 56k, so the bar was drawing the same fact twice
        b.text(v - 700, i, f"{100*v/tot:.0f}%", va="center", ha="right",
               fontsize=6.4, color="white", fontweight="bold")
        b.text(-900, i, SHORT.get(lab, lab), va="center", ha="right",
               fontsize=7.0, color=INK)
    b.set_yticks([]); b.set_xlim(0, tot * 1.02); b.set_ylim(-0.62, len(rows) - 0.38)
    b.set_xticks([0, 20000, 40000, 56012]); b.set_xticklabels(["0", "20k", "40k", "56k"])
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

def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT) as pdf:
        for fn in (figure1, figure2, figure3):
            fig = fn(); pdf.savefig(fig); plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
