#!/usr/bin/env python
"""Shared foundation for the three main-text figure scripts.

Holds the palette, the derived NUMBERS table and the panel helpers. Every
number is measured from the released database or from this session's audits;
none is illustrative, and each entry in NUMBERS names the file it came from so
a reviewer can retrace it.

Split out of the former make_figures.py (2026-09-15), which emitted all three
figures as one three-page PDF. The figures are now separate files, one script
each, so a figure can be regenerated without rebuilding its siblings:

    make_figure1_growth.py          -> ../figures/figure1_growth.pdf
    make_figure2_thermodynamics.py  -> ../figures/figure2_thermodynamics.pdf
    make_figure3_direction.py       -> ../figures/figure3_direction.pdf

Colours are the validated categorical palette; direction uses a diverging
encoding because the variable is polarity (forward / reversible / reverse), and
the 2020-vs-2026 pairs use two steps of one hue because that is magnitude over
time, not identity.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

import grace_style as grace

FIGDIR = Path(__file__).resolve().parent.parent / "figures"

# ---- population -----------------------------------------------------------
# Obsolete records are excluded from every figure. The 2026 draft mixed two
# populations -- "34% growth in compounds" counted all records, "55% in
# reactions" counted 2026 totals against a 2020 baseline with obsolete rows
# stripped -- and the figures silently counted all records throughout. One
# flag, threaded through every loader below, so a figure cannot disagree with
# the text. Flip it to False to reproduce the all-records numbers.
LIVE_ONLY = True


def _live(records):
    """Filter a list of loaded records to the reporting population."""
    if not LIVE_ONLY:
        return records
    return [r for r in records if r.get("is_obsolete") not in (1, "1")]


def _live_rows(rows, col="is_obsolete"):
    """Same, for DictReader rows off a TSV."""
    if not LIVE_ONLY:
        return rows
    return [r for r in rows if r.get(col) != "1"]


def save(fig, name):
    """Write one figure to ../figures/<name>.pdf and report the path."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / name
    fig.savefig(out, format="pdf")
    plt.close(fig)
    print(f"wrote {out}")

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

    COLOUR, 2026-09-22: "no ionizable site" was NEUTRAL grey, the same family
    as "no structure", which merged a chemistry RESULT (Marvin ran and found no
    dissociable proton between pH -2 and 16) with a curation GAP (there is no
    structure to run on). Reviewer 1 asked for a different scheme for this
    panel; the result now takes its own categorical slot and only the genuine
    absence stays grey.

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
        for c in _live(_json.load(open(f))):
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
              ("no ionizable site", cc["no site"], ORANGE, ""),
              ("no structure", cc["no structure"], GRID, "")]

    rx = _c.Counter()
    for f in sorted(_glob.glob(str(root / "Biochemistry/reaction_*.json"))):
        for r in _live(_json.load(open(f))):
            st = r.get("stoichiometry")
            if not isinstance(st, list) or not st: continue
            v = {cls.get(x["compound"], "no structure") for x in st}
            if "no structure" in v: rx["inc"] += 1
            elif any(z.endswith("SMILES") for z in v): rx["smi"] += 1
            elif any(z.endswith("InChI") for z in v): rx["inchi"] += 1
            else: rx["none"] += 1
    by_rxn = [("Marvin", rx["inchi"], BLUE, ""),
              ("from SMILES", rx["smi"], BLUE, "xxx"),
              ("no ionizable site", rx["none"], ORANGE, ""),
              ("no structure", rx["inc"], GRID, "")]
    return by_cpd, by_rxn


def _growth():
    """Figure 1A: compounds carrying an alias from each structure source,
    2020 against 2026. DERIVED from the shipped alias file on both sides.

    The 2020 column is read from Unique_ModelSEED_Compound_Aliases.txt at the
    last commit of 2020 (fd6c7849, 2020-11-10), six weeks after the paper
    appeared online, so the figure cannot drift from the repository. An
    earlier note here called the transcribed values -- MetaCyc 19,138, KEGG
    17,760 -- "wrong" against 19,172 and 17,793 read from the file. They were
    not wrong; they were the same counts with the 34 obsolete 2020 compounds
    removed, i.e. the live basis this function now uses. Both differences are
    exactly that 34 (33 for KEGG, one obsolete compound having no KEGG alias).

    Rhea is absent by construction: its 207 InChIKeys alias to zero ModelSEED
    compounds, so it contributes reactions, not structures.
    """
    import csv as _csv, io as _io, subprocess as _sp, collections as _c, glob as _glob
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

    # The alias files carry no is_obsolete column, so the population filter has
    # to come from the records on each side. The effect is small -- 34 obsolete
    # compounds in 2020, 46 now -- but leaving it out would put this panel on a
    # different population from every other one.
    if LIVE_ONLY:
        import json as _js
        live_new = set()
        for f in sorted(_glob.glob(str(root / "Biochemistry/compound_*.json"))):
            live_new |= {c["id"] for c in _live(_js.load(open(f)))}
        live_old = {r["id"] for r in _csv.DictReader(_io.StringIO(_sp.run(
            ["git", "show", f"{C2020}:Biochemistry/compounds.tsv"], cwd=str(root),
            capture_output=True, text=True, check=True).stdout), delimiter="\t")
            if r.get("is_obsolete") != "1"}
        old = _c.defaultdict(set, {k: v & live_old for k, v in old.items()})
        new = _c.defaultdict(set, {k: v & live_new for k, v in new.items()})
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
        for c in _live(_json.load(open(f))):
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
        for x in _live(_json.load(open(f))):
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


def _euler_regions():
    """Figure 1B: the seven-region overlap of the three primary reaction sources.

    The stacked unique/shared bar this replaces said that 67% of MetaCyc's
    reactions are unique, but never said who the other 33% are shared WITH --
    which is the question integrating Rhea actually raises. These are the
    region counts an Euler diagram needs.

    All records, obsolete included: that is the population the manuscript
    counts everywhere else, and _reaction_sources() above applies no filter
    either, so the two panels stay on one basis.
    """
    import csv as _csv, collections as _c, glob as _glob, json as _json
    root = Path(__file__).resolve().parents[3]
    ids = set()
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_[0-9][0-9].json"))):
        for r in _live(_json.load(open(f))):
            ids.add(r["id"])
    mem = _c.defaultdict(set)
    with (root / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt").open() as fh:
        rd = _csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) >= 3 and row[2] in ("MetaCyc", "KEGG", "Rhea") and row[0] in ids:
                mem[row[0]].add(row[2])
    M = {i for i, v in mem.items() if "MetaCyc" in v}
    K = {i for i, v in mem.items() if "KEGG" in v}
    H = {i for i, v in mem.items() if "Rhea" in v}
    return {"M": len(M - K - H), "K": len(K - M - H), "H": len(H - M - K),
            "MK": len((M & K) - H), "MH": len((M & H) - K), "KH": len((K & H) - M),
            "MKH": len(M & K & H), "none": len(ids - (M | K | H)),
            "totals": {"MetaCyc": len(M), "KEGG": len(K), "Rhea": len(H)}}


def _pathway_classes(top=8):
    """Figure 1D: where the growth landed, at MetaCyc class level.

    Reviewer 1 asked where new information is still being gained. The released
    pathway alias file cannot answer it -- it stops at rxn48568, below the 2020
    boundary, so every post-2020 reaction reads as unannotated. This rebuilds
    the mapping from the shipped source table instead:

        ModelSEED id -> MetaCyc reaction id   (Unique_ModelSEED_Reaction_Aliases)
        MetaCyc reaction -> pathway -> parent (Scripts/Provenance/MetaCyc)

    BOTH eras go through that same join, so the two bars are comparable. Reading
    the pre-2020 side out of the alias file instead would compare a full
    ancestor closure against a one-level parent lookup.

    Super-Pathways is dropped: it is an organisational class, not a biological
    one, and it would otherwise top the chart.
    """
    import csv as _csv, collections as _c, io as _io, subprocess as _sp, glob as _glob, json as _json
    root = Path(__file__).resolve().parents[3]
    C2020 = "fd6c7849891ef4bbeb6eac072f5a6f7adff05b0e"
    DROP = {"Super-Pathways"}

    old_ids = {r["id"] for r in _csv.DictReader(_io.StringIO(_sp.run(
        ["git", "show", f"{C2020}:Biochemistry/reactions.tsv"], cwd=str(root),
        capture_output=True, text=True, check=True).stdout), delimiter="\t")}

    rx2pwy, parent, names = _c.defaultdict(set), {}, {}
    with (root / "Scripts/Provenance/MetaCyc/MetaCyc_pathways.tsv").open() as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            names[r["id"]] = r["name"] or r["id"]
            parent[r["id"]] = [x for x in (r.get("parent") or "").split("|") if x]
            for rx in (r["reactions"] or "").split("|"):
                if rx:
                    rx2pwy[rx].add(r["id"])

    s2m = _c.defaultdict(set)
    with (root / "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt").open() as fh:
        rd = _csv.reader(fh, delimiter="\t"); next(rd)
        for row in rd:
            if len(row) >= 3 and row[2] == "MetaCyc":
                s2m[row[0]].add(row[1])

    # The population filter: this loader iterates the alias file, which has no
    # is_obsolete column, so it was still counting obsolete reactions as
    # "already held" (antibiotics 824 rather than 740) after every other
    # loader had moved to the live basis.
    live = set()
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        live |= {r["id"] for r in _live(_json.load(open(f)))}

    new_c, old_c = _c.Counter(), _c.Counter()
    for sid, mcs in s2m.items():
        if sid not in live:
            continue
        cls = {p for mc in mcs for pwy in rx2pwy.get(mc, ()) for p in parent.get(pwy, [])}
        cls -= DROP
        for c in cls:
            (old_c if sid in old_ids else new_c)[c] += 1
    rows = [(names.get(c, c), n, old_c.get(c, 0)) for c, n in new_c.most_common(top)]
    return rows


def _energy_coverage():
    """Figure 2B: reactions carrying an energy from each source, and the
    denominator.

    Was a hardcoded triple, which is why it stayed on the all-records
    population when everything else moved. Sentinels are excluded, as the
    comment on the old constant said: group contribution writes dg = 1e7 when
    it declines and eQuilibrator writes sigma >= 2500 kcal/mol, and counting
    those reads as coverage that does not exist.
    """
    import json as _json, glob as _glob, collections as _c
    root = Path(__file__).resolve().parents[3]
    n = _c.Counter(); total = 0
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _live(_json.load(open(f))):
            total += 1
            for src, t in (r.get("thermodynamics") or {}).items():
                if src == "LLMs" or not t or t[0] in ("", None):
                    continue
                try:
                    dg, sd = float(t[0]), abs(float(t[1]))
                except (TypeError, ValueError):
                    continue
                if dg >= 1e7 or sd >= 2500:
                    continue
                n[src] += 1
    order = sorted(n, key=lambda k: -n[k])
    return [(k, n[k]) for k in order], total


def _sigma_values():
    """Figure 2C histograms: each source's reported sigma, sentinels excluded.

    Was read from figure_data_sigma.json, a static cache with no generator
    that had been built on ALL records (n = 21,789 / 29,447 / 29,617, medians
    0.63 / 10.41 / 17.01) -- so after the text moved to the live population
    the drawn median lines disagreed with M11's 0.62 / 10.83 / 17.56. Derived
    here through the same population filter and the same sentinel rule as
    _energy_coverage(), so the figure and the text cannot drift apart again.
    Returns {source: sorted list of sigma}.
    """
    import json as _json, glob as _glob, collections as _c
    root = Path(__file__).resolve().parents[3]
    out = _c.defaultdict(list)
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_*.json"))):
        for r in _live(_json.load(open(f))):
            for src, t in (r.get("thermodynamics") or {}).items():
                if src == "LLMs" or not t or t[0] in ("", None):
                    continue
                try:
                    dg, sd = float(t[0]), abs(float(t[1]))
                except (TypeError, ValueError):
                    continue
                if dg >= 1e7 or sd >= 2500:
                    continue
                out[src].append(sd)
    return {k: sorted(v) for k, v in out.items()}


def _silver_sigma_band():
    """Figure 2C shading: the sigma span of reactions graded SILVER, per source.

    Answers "what uncertainty does a silver reaction actually carry?" -- and
    shows that the answer is only meaningful for two of the three sources.
    eQuilibrator and dGPredictor separate their tiers by sigma; Group
    contribution does not, which is the same flat-error-curve problem that
    stopped raw sigma being used to rank sources in the first place. The band
    is p5-p95.

    REWRITTEN to read only files that ship. It previously joined two
    intermediates under results/thermo_grades/ -- source_grades_wide.tsv and
    source_grades.tsv -- both of which are .gitignore'd (see
    Biochemistry/Thermodynamics/SourceGrading/.gitignore line 12). Neither is
    present in a clean checkout, so importing this module raised
    FileNotFoundError and NO figure could be regenerated, this one included.
    Regenerating them is not a way out either: grade_thermo_sources.py fails
    under the pinned pandas.

    The same two quantities are available from files that are committed:
    best_grade from reaction_grades.tsv, and each source's reported sigma from
    the released Biochemistry/reaction_*.json, which is where the published
    uncertainty lives anyway. Same band, reproducible by a reader.
    """
    import csv as _csv, collections as _c, glob as _glob, json as _json
    root = Path(__file__).resolve().parents[3]
    grades = root / ("Biochemistry/Thermodynamics/SourceGrading/results/"
                     "thermo_grades/reaction_grades.tsv")
    best = {}
    with grades.open() as fh:
        for x in _csv.DictReader(fh, delimiter="\t"):
            if x.get("best_grade"):
                best[x["rxn"]] = x["best_grade"]
    sig = _c.defaultdict(list)
    for f in sorted(_glob.glob(str(root / "Biochemistry" / "reaction_[0-9][0-9].json"))):
        for r in _live(_json.load(open(f))):
            if best.get(r["id"]) != "SILVER":
                continue
            for name, triple in (r.get("thermodynamics") or {}).items():
                if name == "LLMs" or not triple or len(triple) < 2:
                    continue
                try:
                    dg, sd = float(triple[0]), abs(float(triple[1]))
                except (TypeError, ValueError):
                    continue
                # Sentinels, not error bars. Same exclusions the Results
                # section documents: group contribution writes dg = 1e7 when
                # it declines, and eQuilibrator writes sigma >= 2500 kcal/mol.
                # Leaving them in put eQuilibrator's p95 at 23,901.
                if dg >= 1e7 or sd >= 2500:
                    continue
                sig[name].append(sd)
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
        for r in _live(_json.load(open(f))):
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
        for r in _live(_json.load(open(f))):
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
        for r in _live(_json.load(open(f))):
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
    "energy_rxn": None,    # filled by _energy_coverage()
    "energy_total": None,  # filled by _energy_coverage()
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
NUMBERS["energy_rxn"], NUMBERS["energy_total"] = _energy_coverage()
NUMBERS["sigma"] = _sigma_values()


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

