#!/usr/bin/env python
"""Write the reaction-level evidence grade into `thermo-evidence`, and the
recommended direction into `reversibility`.

WHAT GOES WHERE
    thermo-evidence = {grade, assessment, cross-source, source}
        A SIBLING of `thermodynamics`, not a member of it: that dict is strictly
        source-keyed and a consumer iterating it must not meet a key that is not
        a source.
        One grade per REACTION, formed from its collective evidence. The three
        labels answer: how good is the evidence (grade), what did the deciding
        source claim about itself (self_assessment), and what did the other
        sources make of it (cross_source, absent when the check resolved
        neither way). `source` names the source the grade came from.

    reversibility = the recommended direction, by precedence
        eQuilibrator > dGPredictor > Group contribution -- the first of those
        with a callable direction wins. GC is last because it is the legacy
        source: its reported error does not separate its own grade tiers and it
        overstates by ~2.2x. It still recommends for 4,706 reactions no modern
        source can call.

THE PREVIOUS `reversibility` WAS THE 2020 CANONICAL FIELD (eQuilibrator primary,
Group Contribution fallback, historical rule set). It is regenerable at any time
with Apply_2020_Reversibility_Policy.py, which is how to undo this.
"""
import argparse, csv, glob, json, os, sys
from collections import Counter, defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GRADES = os.path.join(REPO, "Biochemistry/Thermodynamics/SourceGrading/results/thermo_grades")
KEY = {"EQ": "eQuilibrator", "GC": "Group contribution",
       "DGP": "dGPredictor", "TECRDB": "TECRDB"}
VERDICTS = ("corroborated", "outvoted", "unpaired")
ABBREV = {"eQuilibrator": "eQ", "Group contribution": "GC",
          "dGPredictor": "dGP", "TECRDB": "TECRDB"}
PRECEDENCE = ["eQuilibrator", "dGPredictor", "Group contribution"]


def load_grades():
    reason = defaultdict(dict)
    with open(os.path.join(GRADES, "source_grades.tsv")) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            reason[x["rxn"]][x["source"]] = x["reason"]
    out = {}
    with open(os.path.join(GRADES, "reaction_grades.tsv")) as fh:
        for x in csv.DictReader(fh, delimiter="\t"):
            if not x["best_grade"]:
                continue
            src = KEY.get(x["best_source"], x["best_source"])
            rsn = reason.get(x["rxn"], {}).get(src, "")
            if rsn == "measured":
                conf, cross = "measured", None
            else:
                conf, cross = rsn, None
                for v in VERDICTS:
                    if rsn.endswith("-" + v):
                        conf, cross = rsn[:-(len(v) + 1)], v
                        break
            ev = {"grade": x["best_grade"].lower(), "assessment": conf,
                  "source": ABBREV.get(src, src)}
            if cross:
                ev["cross-source"] = cross
            out[x["rxn"]] = ev
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="report only, write nothing")
    a = ap.parse_args()

    ev = load_grades()
    stats = Counter()
    example = None
    for path in sorted(glob.glob(os.path.join(REPO, "Biochemistry", "reaction_*.json"))):
        data = json.load(open(path))
        for r in data:
            th = r.get("thermodynamics") or {}
            ops = {s: v[2] for s, v in th.items()
                   if isinstance(v, list) and len(v) > 2 and v[2] in ("<", ">", "=")}
            rec = next((ops[s] for s in PRECEDENCE if s in ops), "?")
            src = next((s for s in PRECEDENCE if s in ops), None)
            e, d = ops.get("eQuilibrator"), ops.get("dGPredictor")
            conflict = bool(e and d and e != d)

            th.pop("evidence", None)          # never inside thermodynamics
            r.pop("direction_conflict", None)  # retired
            if r["id"] in ev:
                r["thermo-evidence"] = ev[r["id"]]
                stats["graded"] += 1
            else:
                r.pop("thermo-evidence", None)
            stats["rec_" + (src or "none")] += 1
            r["reversibility"] = rec
            if example is None and r["id"] in ev and "cross-source" in ev[r["id"]]:
                example = json.dumps(
                    {k: r[k] for k in ("id", "name", "reversibility",
                                       "thermodynamics", "thermo-evidence")},
                    indent=2)[:1000]
        if not a.dry_run:
            with open(path, "w") as fh:
                json.dump(data, fh, indent=4, sort_keys=True)
    print(f"  reactions carrying an evidence grade : {stats['graded']:,}")
    for s in PRECEDENCE + ["none"]:
        print(f"  recommendation from {s:22s}: {stats['rec_' + s]:,}")
    print("\n  EXAMPLE:")
    print("\n".join("    " + l for l in (example or "").splitlines()))
    if a.dry_run:
        print("\n  DRY RUN -- nothing written.")


if __name__ == "__main__":
    main()
