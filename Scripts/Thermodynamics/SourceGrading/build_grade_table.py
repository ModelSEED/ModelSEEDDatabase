#!/usr/bin/env python3
"""Emit the LaTeX body of the grades table (Supplementary Table S1).

Written 2026-09-14 because no generator existed: the table in the supplement had
been assembled by hand and had drifted from the data (~1,550 reactions a tier
out, and a Conflict column that could not be reproduced under any definition).
Regenerate with this and paste, or diff against the committed table.

Rows are (best_grade, self-assessment, cross-source), per reaction, taking the
reason of the source that earned the reaction's best grade.

Two collapses, per Sam 2026-09-14 -- a verdict is named only where it changed
the grade:
  * corroborated folds into '---' for self-certain / self-confident
    (p_ok alone already set the tier; corroboration lifts bronze only)
  * unpaired folds into '---' for self-confident / unconfident
    (grade_thermo_sources.py: `g[is_unpaired & (g == GOLD)] = SILVER`,
     so unpaired demotes gold and nothing else)

CONFLICT counts reactions where eQuilibrator and dGPredictor both commit to an
IRREVERSIBLE direction and those directions are opposed ('>' vs '<'). A source
saying '>' against '=' is a difference of confidence, not a contradiction.
"""
import csv, json, glob, collections, os
from pathlib import Path

ROOT = Path(os.environ.get("MSDB_ROOT", Path(__file__).resolve().parents[3]))
G = ROOT / "Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades"
FULL = {"EQ": "eQuilibrator", "DGP": "dGPredictor",
        "GC": "Group contribution", "openTECR": "openTECR"}
SELF = ["measured", "self-certain", "self-confident", "unconfident"]
SORT_SELF = {s: i for i, s in enumerate(SELF)}
SORT_CROSS = {"corroborated": 0, "disputed": 1, "unpaired": 2, "---": 3}
SORT_GRADE = {"gold": 0, "silver": 1, "bronze": 2}
# Notation, fixed 2026-09-15. Every condition is its own math group and the
# separator is always "; ", outside math. Previously the p_ok clause and the
# R/z clause were joined by ";" while R and z were joined by an in-math ","
# -- two separators for one job -- and \! negative thin spaces jammed each
# relation against its operands, so "$R\!>\!2, z\!>\!2$" set as a cramped
# "R>2,z>2". Keep conditions atomic; do not reintroduce \!.
PARAM = {
    ("gold", "measured"):       r"stereo-exact anchor",
    ("gold", "self-certain"):   r"$p_{\text{ok}} \ge 0.90$",
    ("silver", "self-certain"): r"$p_{\text{ok}} \ge 0.90$",
    ("silver", "self-confident"): r"$0.70 \le p_{\text{ok}} < 0.90$",
    ("silver", "unconfident"):  r"$p_{\text{ok}} < 0.70$",
    ("bronze", "self-confident"): r"$0.70 \le p_{\text{ok}} < 0.90$",
    ("bronze", "unconfident"):  r"$p_{\text{ok}} < 0.70$",
}
EXTRA = {"corroborated": r"; $R \le 2$; $z \le 2$",
         "disputed":     r"; $R > 2$; $z > 2$",
         "unpaired":     r"; $n_{\text{src}} = 1$", "---": ""}


def split_reason(rs):
    for s in SELF:
        if rs == s:
            return s, "---"
        if rs.startswith(s + "-"):
            return s, rs[len(s) + 1:]
    raise ValueError(f"unparsable reason {rs!r}")


def main():
    rg = {r["rxn"]: r for r in csv.DictReader(open(G / "reaction_grades.tsv"), delimiter="\t")}
    sg = collections.defaultdict(dict)
    for r in csv.DictReader(open(G / "source_grades.tsv"), delimiter="\t"):
        sg[r["rxn"]][r["source"]] = r

    direction = {}
    for f in sorted(glob.glob(str(ROOT / "Biochemistry/reaction_*.json"))):
        for x in json.load(open(f)):
            t = x.get("thermodynamics") or {}
            direction[x["id"]] = (t.get("eQuilibrator", ["", "", ""])[2],
                                  t.get("dGPredictor", ["", "", ""])[2])

    n = collections.Counter()
    conflict = collections.Counter()
    for rx, r in rg.items():
        grade = r["best_grade"]
        if not grade:
            continue
        row = sg[rx].get(FULL.get(r["best_source"], r["best_source"]))
        if row is None:
            raise SystemExit(f"{rx}: best_source {r['best_source']} absent from source_grades.tsv")
        s, x = split_reason(row["reason"])
        if x == "corroborated" and s in ("self-certain", "self-confident"):
            x = "---"
        if x == "unpaired" and s in ("self-confident", "unconfident"):
            x = "---"
        key = (grade.lower(), s, x)
        n[key] += 1
        e, d = direction.get(rx, ("", ""))
        if e in (">", "<") and d in (">", "<") and e != d:
            conflict[key] += 1

    keys = sorted(n, key=lambda k: (SORT_GRADE[k[0]], SORT_SELF[k[1]], SORT_CROSS[k[2]]))
    fmt = lambda v: f"{v:,}".replace(",", "{,}")
    prev = None
    for k in keys:
        if prev is not None and k[0] != prev:
            print(r"\midrule")
        param = PARAM[(k[0], k[1])] + EXTRA[k[2]]
        print(f"{k[0]:<6} & {k[1]:<14} & {k[2]:<12} & {param} & {fmt(n[k])} & {fmt(conflict[k])} \\\\")
        prev = k[0]
    print(r"\midrule")
    print(r"\multicolumn{4}{l}{\textbf{Total}} & \textbf{%s} & \textbf{%s} \\"
          % (fmt(sum(n.values())), fmt(sum(conflict.values()))))


if __name__ == "__main__":
    main()
