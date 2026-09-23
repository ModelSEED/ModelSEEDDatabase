#!/usr/bin/env python
"""Every count the manuscript quotes, on both populations, side by side.

The 2026 draft mixed two populations: "34% growth in compounds" counted all
records at both time points, while "55% in reactions" set 2026's totals against
a 2020 baseline with its obsolete rows removed. Neither number was wrong on its
own terms; together they were not comparable.

This prints each quantity on both bases so the conversion can be audited rather
than trusted. ALL RECORDS is what the submitted draft used throughout apart
from that one figure; LIVE drops `is_obsolete` rows.

The definitions matter more than the filter, and they are the ones that
reproduce the draft exactly on the all-records column:

  covers      carries an energy, EXCLUDING each source's refusal sentinel --
              group contribution writes dg = 1e7 when it declines, and
              eQuilibrator writes sigma >= 2500 kcal/mol. Counting the
              sentinels puts eQuilibrator's coverage at 25,175 against the
              21,789 the draft quotes.
  resolves    states a direction OR an explicit reversibility: '>', '<' or '='.
              Restricting this to '>' and '<' turns the draft's 13,289 into
              6,267, which is how the definition announced itself.

Run:  ~/Documents/py_venv/bin/python population_basis_table.py
"""
import csv
import glob
import io
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
C2020 = "fd6c7849891ef4bbeb6eac072f5a6f7adff05b0e"
PREDICTORS = ("eQuilibrator", "dGPredictor", "Group contribution")
STATED = (">", "<", "=")
COMMITTED = (">", "<")


def _load(pattern):
    out = []
    for f in sorted(glob.glob(str(ROOT / "Biochemistry" / pattern))):
        out += json.load(open(f))
    return out


def direction(rxn, source):
    t = (rxn.get("thermodynamics") or {}).get(source)
    return (t[2] or None) if t and len(t) > 2 else None


def energy(rxn, source):
    """The reported energy, or None when the source declined. Sentinels are
    refusal markers, not values."""
    t = (rxn.get("thermodynamics") or {}).get(source)
    if not t or t[0] in ("", None):
        return None
    try:
        dg, sd = float(t[0]), abs(float(t[1]))
    except (TypeError, ValueError):
        return None
    if dg >= 1e7 or sd >= 2500:
        return None
    return dg


def row(label, a, b, paper=""):
    fa = f"{a:,}" if isinstance(a, int) else f"{a:.1f}%"
    fb = f"{b:,}" if isinstance(b, int) else f"{b:.1f}%"
    print(f"  {label:<44}{fa:>12}{fb:>12}   {paper}")


def main():
    rxns = _load("reaction_[0-9][0-9].json")
    cpds = _load("compound_[0-9][0-9].json")
    live_r = [r for r in rxns if r.get("is_obsolete") != 1]
    live_c = [c for c in cpds if c.get("is_obsolete") != 1]

    old = subprocess.run(["git", "show", f"{C2020}:Biochemistry/reactions.tsv"],
                         cwd=str(ROOT), capture_output=True, text=True,
                         check=True).stdout
    old_r = list(csv.DictReader(io.StringIO(old), delimiter="\t"))
    oldc = subprocess.run(["git", "show", f"{C2020}:Biochemistry/compounds.tsv"],
                          cwd=str(ROOT), capture_output=True, text=True,
                          check=True).stdout
    old_c = list(csv.DictReader(io.StringIO(oldc), delimiter="\t"))
    old_r_live = [r for r in old_r if r.get("is_obsolete") != "1"]
    old_c_live = [c for c in old_c if c.get("is_obsolete") != "1"]

    print("=" * 86)
    print(f"  {'quantity':<44}{'ALL RECORDS':>12}{'LIVE':>12}   draft says")
    print("=" * 86)

    print("\nTOTALS  (abstract, M09)")
    row("compounds", len(cpds), len(live_c), "~46,000")
    row("reactions", len(rxns), len(live_r), "~56,000")
    st_a = sum(1 for c in cpds if c.get("smiles") or c.get("inchikey"))
    st_l = sum(1 for c in live_c if c.get("smiles") or c.get("inchikey"))
    row("compounds with a structure", st_a, st_l, "36,943 / ~37,000")
    row("2020 compounds", len(old_c), len(old_c_live), "33,992")
    row("2020 reactions", len(old_r), len(old_r_live), "43,774")
    print(f"  {'GROWTH, compounds':<44}"
          f"{100*(len(cpds)/len(old_c)-1):>11.1f}%{100*(len(live_c)/len(old_c_live)-1):>11.1f}%   34%")
    print(f"  {'GROWTH, reactions':<44}"
          f"{100*(len(rxns)/len(old_r)-1):>11.1f}%{100*(len(live_r)/len(old_r_live)-1):>11.1f}%"
          f"   55% (mixed basis)")

    print("\nTHERMODYNAMIC COVERAGE  (M11 Fig 2B, M12)")
    for s in PREDICTORS:
        a = sum(1 for r in rxns if energy(r, s) is not None)
        b = sum(1 for r in live_r if energy(r, s) is not None)
        row(f"{s} covers", a, b)
    print()
    for s in PREDICTORS:
        ca = [r for r in rxns if energy(r, s) is not None]
        cb = [r for r in live_r if energy(r, s) is not None]
        pa = 100 * sum(1 for r in ca if direction(r, s) in STATED) / len(ca)
        pb = 100 * sum(1 for r in cb if direction(r, s) in STATED) / len(cb)
        row(f"{s} resolves", pa, pb)

    print("\nDIRECTION  (M12)")
    a = sum(1 for r in rxns if any(direction(r, s) in STATED for s in PREDICTORS))
    b = sum(1 for r in live_r if any(direction(r, s) in STATED for s in PREDICTORS))
    row("reactions directed by >=1 source", a, b, "30,157")
    row("reactions directed by none", len(rxns) - a, len(live_r) - b, "25,855")
    for name, pool, tot in (("all", rxns, None), ("live", live_r, None)):
        pass
    eq_a = {r["id"] for r in rxns if direction(r, "eQuilibrator") in STATED}
    dg_a = {r["id"] for r in rxns if direction(r, "dGPredictor") in STATED}
    eq_b = {r["id"] for r in live_r if direction(r, "eQuilibrator") in STATED}
    dg_b = {r["id"] for r in live_r if direction(r, "dGPredictor") in STATED}
    row("dGPredictor resolves, eQuilibrator not", len(dg_a - eq_a), len(dg_b - eq_b), "4,248")
    row("eQuilibrator resolves, dGPredictor not", len(eq_a - dg_a), len(eq_b - dg_b), "13,289")

    def agree(pool):
        e = {r["id"]: direction(r, "eQuilibrator") for r in pool
             if direction(r, "eQuilibrator") in COMMITTED}
        d = {r["id"]: direction(r, "dGPredictor") for r in pool
             if direction(r, "dGPredictor") in COMMITTED}
        both = set(e) & set(d)
        return 100 * sum(1 for i in both if e[i] == d[i]) / len(both)
    row("eQ/dGP agreement where both commit", agree(rxns), agree(live_r), "95.0%")

    print("\nRECOMMENDED DIRECTION  (M12)")
    for src, paper in (("eQ", "21,218"), ("dGP", "4,248"), ("GC", "4,691")):
        a = sum(1 for r in rxns if (r.get("thermo-evidence") or {}).get("source") == src)
        b = sum(1 for r in live_r if (r.get("thermo-evidence") or {}).get("source") == src)
        row(f"recommended from {src}", a, b, paper)

    print("\nEVIDENCE GRADES  (M12, Supplementary Table S1)")
    ga = Counter((r.get("thermo-evidence") or {}).get("grade") for r in rxns)
    gb = Counter((r.get("thermo-evidence") or {}).get("grade") for r in live_r)
    row("graded reactions", sum(v for k, v in ga.items() if k),
        sum(v for k, v in gb.items() if k), "33,099")
    for t, paper in (("gold", "3,434"), ("silver", "18,388"), ("bronze", "11,277")):
        row(f"  {t}", ga[t], gb[t], paper)
    aa = Counter((r.get("thermo-evidence") or {}).get("assessment") for r in rxns)
    ab = Counter((r.get("thermo-evidence") or {}).get("assessment") for r in live_r)
    for t, paper in (("measured", "806"), ("self-certain", "2,628"),
                     ("self-confident", "9,477+161"), ("unconfident", "")):
        row(f"  assessed {t}", aa[t], ab[t], paper)

    print("\n  NOTE on the 806 -> 365 anchors: the 441 records dropped are NOT lost")
    print("  measurements. Every one is flagged obsolete AND carries a linked")
    print("  replacement -- they are duplicate records of the same chemistry, so a")
    print("  single measurement was being counted several times. 365 is the number")
    print("  of DISTINCT anchored reactions, and is the more defensible figure.")

    print("\nATOM MAPPING  (M13)")
    am_a = [r for r in rxns if r.get("atom_mapping")]
    am_b = [r for r in live_r if r.get("atom_mapping")]
    row("reactions with a mapping", len(am_a), len(am_b), "32,877")
    for conf in ("clean", "salvaged"):
        row(f"  {conf}",
            sum(1 for r in am_a if (r["atom_mapping"] or {}).get("confidence") == conf),
            sum(1 for r in am_b if (r["atom_mapping"] or {}).get("confidence") == conf),
            {"clean": "25,058", "salvaged": "7,819"}[conf])
    print(f"  {'  as % of reactions':<44}"
          f"{100*len(am_a)/len(rxns):>11.1f}%{100*len(am_b)/len(live_r):>11.1f}%   59%")

    print("\nTRANSPORT  (M05)")
    ta = [r for r in rxns if r.get("is_transport") == 1]
    tb = [r for r in live_r if r.get("is_transport") == 1]
    row("transport reactions", len(ta), len(tb))
    for pool, lbl in ((ta, "all"), (tb, "live")):
        gold = [r for r in pool if (r.get("thermo-evidence") or {}).get("grade") == "gold"]
        atp = sum(1 for r in gold
                  if {"cpd00002", "cpd00008"} <=
                  {p["compound"] for p in r.get("stoichiometry") or []})
        print(f"    {lbl:<6} gold {len(gold):>6,}   ATP-coupled {100*atp/len(gold):5.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
