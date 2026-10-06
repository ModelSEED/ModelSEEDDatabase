#!/usr/bin/env python
"""#48 -- complete macroscopic pKa ladders for the head of the #31 ranking.

Emits two files.

``ladder_requirements.tsv``
    For each head compound, the ladder SHAPE the magnesium guard in
    ``build_modelseed_cache.py`` demands. The guard is::

        major    = sum(1 for v in ladder if v > 7.0)
        reachable = {i - major + n_h for i in range(len(ladder) + 1)}
        proceed if set(mg_protons) <= reachable

    so with ``lo = min(mg_protons)`` and ``hi = max(mg_protons)`` a ladder is
    admissible iff it carries at least ``n_above = max(0, n_h - lo)`` values
    above pH 7 and at least ``n_below = max(0, hi - n_h)`` values at or below
    it. Both numbers are reported, plus whether the sources we hold satisfy
    them.

``literature_ladders.tsv``
    The curated table. **Every row carries a resolvable citation.** Rows are
    emitted only from sources whose provenance is machine-checkable:

    * ``alberty`` -- BasicBiochemData3, doi:10.1002/0471332607, ladder derived
      by Alberty's own calcpK relation from the tabulated species dfG values.
    * ``iupac``   -- the IUPAC Digitized pKa Dataset (Zenodo 10.5281/
      zenodo.15375522), whose per-row ``ref`` code is expanded here to the full
      Serjeant/Perrin bibliographic entry via ``reference_code_translation.csv``.

    Nothing is written from recollection. A dissociation with no source in
    hand is left out of this file and recorded as a gap in the requirements
    file instead, because a table mixing cited and remembered values is worse
    than a short one.
"""
import csv, io, pickle, sqlite3, sys, zipfile
from pathlib import Path

from paths import (EQ, PKA_DIR, DATA_OUT, THERMO_CACHE, RANKED,
                   IUPAC_ZIP, require)

OUT = DATA_OUT
ZIP = IUPAC_ZIP
CACHE = THERMO_CACHE
RANK = RANKED
TOP_N = 25

sys.path.insert(0, str(EQ / "tools"))
from modelseed_pkas import load_alberty_pkas, cache_seed_identifiers  # noqa: E402

ALBERTY_CITATION = ("R. A. Alberty, Thermodynamics of Biochemical Reactions, "
                    "Wiley, 2003, doi:10.1002/0471332607 (BasicBiochemData3 "
                    "package; ladder derived from tabulated species dfG by "
                    "pKa = [dfG(deprot) - dfG(prot)] / RT ln10, I = 0)")


def iupac_rows_with_citations():
    """seed_id -> list of dicts, each an IUPAC dissociation with a full citation."""
    if not ZIP.exists():
        return {}
    z = zipfile.ZipFile(ZIP)
    refname = next(n for n in z.namelist() if n.endswith("reference_code_translation.csv"))
    refmap = {}
    for r in csv.DictReader(io.StringIO(z.read(refname).decode("utf-8", "ignore"))):
        refmap[r["code"]] = r
    # the per-compound file the pipeline already consumes, which carries the
    # median value, the ref codes and the remark text
    tbl = PKA_DIR / "iupac_v2_3b.tsv"
    if not tbl.exists():
        return {}
    # source book per (seed compound, pka_type) is needed to pick the right
    # column of the translation table; recover it from the raw dataset
    raw = next(n for n in z.namelist() if n.endswith("iupac_high-confidence_v2_3.csv"))
    bysource = {}
    for r in csv.DictReader(io.StringIO(z.read(raw).decode("utf-8", "ignore"))):
        bysource.setdefault((r["pka_type"], r["ref"]), set()).add(r["source"])
    out = {}
    for r in csv.DictReader(tbl.open(), delimiter="\t"):
        cites = []
        for code in [c.strip() for c in r["refs"].split(";") if c.strip()]:
            code = code.split(",")[0].strip()
            ent = refmap.get(code)
            if not ent:
                continue
            srcs = bysource.get((r["pka_type"], code)) or set()
            col = ("serjeant_source" if "serjeant" in srcs else
                   "perrin_source" if "perrin" in srcs else
                   "perrin_supp_source" if "perrin_supp" in srcs else None)
            text = (ent.get(col) or "").strip() if col else ""
            if not text:                      # fall back to whichever is populated
                text = next((ent[c].strip() for c in
                             ("serjeant_source", "perrin_source", "perrin_supp_source")
                             if ent.get(c, "").strip()), "")
            if text:
                cites.append(f"[{code}] {text}")
        if not cites:
            continue
        out.setdefault(r["seed_id"], []).append({
            "step": r["pka_type"], "value": float(r["pka_value"]),
            "temperature_C": r["temperature_C"], "assessment": r["assessment"],
            "remarks": r["remarks"], "n_measurements": r["n_measurements"],
            "spread": r["spread"], "citation": " | ".join(cites),
        })
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    alberty = load_alberty_pkas()
    iupac = iupac_rows_with_citations()
    sm = cache_seed_identifiers(str(CACHE))
    con = sqlite3.connect(str(CACHE))
    rank = list(csv.DictReader(RANK.open(), delimiter="\t"))[:TOP_N]

    req = csv.writer((OUT / "ladder_requirements.tsv").open("w", newline=""), delimiter="\t")
    req.writerow(["rank", "seed_id", "name", "formula", "n_reactions", "n_h",
                  "mg_protons", "need_above_pH7", "need_at_or_below_pH7",
                  "need_total_steps", "have_alberty", "have_iupac",
                  "alberty_guard", "iupac_guard", "resolver_pick",
                  "resolver_guard", "still_missing_steps"])
    HDR = ["seed_id", "name", "step", "pka_value", "kind", "source",
           "temperature_C", "ionic_strength", "assessment",
           "n_measurements", "spread", "citation", "remarks"]
    # Two outputs. The combined table carries IUPAC values and is GITIGNORED:
    # the IUPAC Digitized pKa Dataset is CC-BY-NC-4.0 and this database is
    # redistributed without a non-commercial restriction, so those values are
    # consumed locally and never committed (same rule as iupac_*.tsv, and the
    # decision recorded in task #32). The open table is Alberty only and is
    # safe to redistribute. Citations are printable in both cases; it is the
    # measured values that are encumbered.
    _all = (OUT / "literature_ladders.tsv").open("w", newline="")
    _open = (OUT / "literature_ladders_open.tsv").open("w", newline="")
    lit_all, lit_open = csv.writer(_all, delimiter="\t"), csv.writer(_open, delimiter="\t")
    lit_all.writerow(HDR); lit_open.writerow(HDR)

    class _Fan:
        def writerow(self, row):
            lit_all.writerow(row)
            if row[5] == "alberty":
                lit_open.writerow(row)
    lit = _Fan()

    for i, r in enumerate(rank, 1):
        sid = r["seed_id"]
        ent = sm.get(sid)
        if not ent:
            continue
        cid = ent[0]
        row = con.execute("select atom_bag from compounds where id=?", (cid,)).fetchone()
        bag = row[0]
        bag = pickle.loads(bag) if isinstance(bag, (bytes, bytearray)) else (bag or {})
        n_h = bag.get("H", 0) if isinstance(bag, dict) else 0
        mgp = sorted({x[0] for x in con.execute(
            "select number_protons from magnesium_dissociation_constant "
            "where compound_id=?", (cid,))})
        if mgp:
            need_above = max(0, n_h - min(mgp))
            need_below = max(0, max(mgp) - n_h)
        else:
            need_above = need_below = 0
        need_total = need_above + need_below

        a = alberty.get(sid, [])
        u = [d["value"] for d in iupac.get(sid, [])]

        def passes(lad):
            if not lad:
                return False
            major = sum(1 for v in lad if v > 7.0)
            reach = {j - major + n_h for j in range(len(lad) + 1)}
            return set(mgp) <= reach if mgp else True

        def verdict(lad):
            if not mgp:
                return "n/a (no Mg data)"
            if not lad:
                return "no ladder"
            return "PASS" if passes(lad) else "FAIL"

        # the resolver takes exactly one source per compound, alberty first
        pick = "alberty" if a else ("iupac" if u else "none")
        best = a or u
        ok = bool(mgp) and passes(best) if best else False
        short = ""
        if mgp and not ok:
            short = (f"need >={need_above} above pH7 and >={need_below} at/below; "
                     f"alberty {sum(1 for v in a if v>7)} above/"
                     f"{len(a)-sum(1 for v in a if v>7)} below, "
                     f"iupac {sum(1 for v in u if v>7)} above/"
                     f"{len(u)-sum(1 for v in u if v>7)} below")
        req.writerow([i, sid, r["name"], r["formula"], r["n_reactions"], n_h,
                      ";".join(map(str, mgp)) or "-", need_above, need_below,
                      need_total, len(a), len(u), verdict(a), verdict(u),
                      pick, verdict(best), short])

        for n, v in enumerate(a, 1):
            lit.writerow([sid, r["name"], f"alberty_step{n}", f"{v:.4f}", "macroscopic",
                          "alberty", "25", "0 (I=0)", "derived", "", "",
                          ALBERTY_CITATION, "Alberty tabulates only species "
                          "significant between roughly pH 5 and 9"])
        for d in iupac.get(sid, []):
            lit.writerow([sid, r["name"], d["step"], f"{d['value']:.4f}", "macroscopic",
                          "iupac", d["temperature_C"], "see remarks",
                          d["assessment"], d["n_measurements"], d["spread"],
                          d["citation"], d["remarks"]])
    _all.close()
    _open.close()
    print("wrote", OUT / "ladder_requirements.tsv")
    print("wrote", OUT / "literature_ladders.tsv", "(gitignored: carries IUPAC values)")
    print("wrote", OUT / "literature_ladders_open.tsv", "(Alberty only, redistributable)")


if __name__ == "__main__":
    main()
