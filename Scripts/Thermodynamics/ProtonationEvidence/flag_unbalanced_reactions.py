#!/usr/bin/env python
"""Flag reactions whose protonation model is unreliable in the shipped build.

Substituting a compound's pKa source changes how many of its sites lie above
the reported pH, and each such site puts a proton on the major microspecies.
Where a reaction's partners are substituted too the change cancels; where they
are not, it does not, and dG'0 moves by RT ln(10) x pH = 9.55 kcal/mol per
unbalanced proton.

This compares the shipped build against a reference build and lists every
reaction whose energy moved more than --threshold, recording whether the
unbalanced-enumeration signature is present.

IT IS A RISK MARKER, NOT A VERDICT. The signature can hold on a reaction that
nonetheless came out right: rxn08713 satisfies it and lands on the
independently-derived Marvin-era value. And rxn08789 moved more than the
threshold WITHOUT the signature and is a correction, not damage. The flag says
"do not read this at face value", not "this is wrong".

Nor is the reference build a safe fallback. On these reactions both partners
carried the same spurious sites there and the errors cancelled, which is not
the same as being correct.
"""
import argparse
import csv
import glob
import json
import pickle
import re
import sqlite3
import sys
from pathlib import Path

from paths import DATA_OUT, EQ, REPO, require

sys.path.insert(0, str(EQ / "tools"))
from modelseed_pkas import cache_seed_identifiers  # noqa: E402

CPD = re.compile(r"cpd\d{5}")


def sites_above(db, ph=7.0):
    seeds = cache_seed_identifiers(str(db))
    con = sqlite3.connect(str(db))
    out = {}
    for seed, (cid, _k) in seeds.items():
        row = con.execute(
            "select dissociation_constants from compounds where id=?", (cid,)).fetchone()
        if not row:
            continue
        v = row[0]
        try:
            v = pickle.loads(v) if isinstance(v, (bytes, bytearray)) else v
        except Exception:
            continue
        out[seed] = sum(1 for x in (v or []) if x > ph)
    con.close()
    return out


def load(path):
    with open(path) as fh:
        return {r["reaction_id"]: r
                for r in csv.DictReader((l for l in fh if not l.startswith("#")),
                                        delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shipped-cache", type=Path,
                    default=EQ / "data/cache_final/compounds.sqlite")
    ap.add_argument("--shipped-energies", type=Path,
                    default=REPO / "Biochemistry/Thermodynamics/eQuilibrator"
                                   "/ModelSEED_Reaction_Energies.tsv")
    ap.add_argument("--reference-cache", type=Path,
                    default=EQ / "data/prev_20260903/cache_final/compounds.sqlite")
    ap.add_argument("--reference-energies", type=Path,
                    default=EQ / "data/prev_20260903/reactions_final.tsv")
    ap.add_argument("--threshold", type=float, default=10.0,
                    help="kcal/mol; a move larger than this is flagged")
    ap.add_argument("--out", type=Path,
                    default=DATA_OUT / "unbalanced_enumeration_reactions.tsv")
    a = ap.parse_args()

    for p, what in ((a.shipped_cache, "shipped cache"),
                    (a.reference_cache, "reference cache"),
                    (a.shipped_energies, "shipped energies"),
                    (a.reference_energies, "reference energies")):
        require(p, what)

    before, after = sites_above(a.reference_cache), sites_above(a.shipped_cache)
    ref, ship = load(a.reference_energies), load(a.shipped_energies)

    names = {}
    for f in glob.glob(str(REPO / "Biochemistry" / "compound_*.json")):
        for c in json.load(open(f)):
            names[c["id"]] = (c.get("name") or "")[:44]

    rows = []
    for rid, r in ship.items():
        b = ref.get(rid)
        if not b or r["status"] != "ok" or b["status"] != "ok":
            continue
        try:
            d = float(r["dg_prime_kcal_per_mol"]) - float(b["dg_prime_kcal_per_mol"])
        except (TypeError, ValueError):
            continue
        if abs(d) <= a.threshold:
            continue
        cs = set(CPD.findall(r.get("formula") or ""))
        moved = sorted(c for c in cs if after.get(c, 0) != before.get(c, 0))
        stay = sorted(c for c in cs if c not in moved and before.get(c, 0) > 0)
        big = max(stay, key=lambda c: before.get(c, 0)) if stay else ""
        rows.append({
            "reaction_id": rid, "name": (r.get("name") or "")[:60],
            "dg_reference_kcal": b["dg_prime_kcal_per_mol"],
            "dg_shipped_kcal": r["dg_prime_kcal_per_mol"],
            "shift_kcal": f"{d:.2f}",
            "net_site_change": sum(after.get(c, 0) - before.get(c, 0) for c in moved),
            "flag": "unbalanced_enumeration" if (moved and stay)
                    else "large_shift_mechanism_unidentified",
            "compounds_resubstituted": ";".join(moved),
            "partners_left_unchanged": ";".join(stay),
            "largest_partner": big, "largest_partner_name": names.get(big, "")})

    rows.sort(key=lambda r: -abs(float(r["shift_kcal"])))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        fh.write(f"""# Reactions whose energy moved more than {a.threshold} kcal/mol between the
# reference and shipped builds. Treat every one as LOW CONFIDENCE.
#
# RISK MARKER, NOT A VERDICT:
#   unbalanced_enumeration  a compound's pKa source was substituted while a
#       partner carrying sites above pH 7 was not, so spurious sites no longer
#       cancel (~9.55 kcal/mol per unbalanced proton). rxn08713 satisfies this
#       and still lands on the independently-derived Marvin-era value, so the
#       signature over-selects.
#   large_shift_mechanism_unidentified  moved without that signature. Includes
#       rxn08789, which moved TOWARD the Marvin-era value and is a correction.
#
# The reference value is NOT a safe fallback: there both partners carried the
# same spurious sites and the errors cancelled, which is not accuracy.
# NONE of these reactions has a TECRDB match, so neither value can be checked
# against experiment.
#
# shipped:   {a.shipped_cache}
# reference: {a.reference_cache}
""")
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    import collections
    print(f"flagged {len(rows)} reactions -> {a.out}")
    for k, v in collections.Counter(r["flag"] for r in rows).most_common():
        print(f"   {k:38s} {v}")


if __name__ == "__main__":
    main()
