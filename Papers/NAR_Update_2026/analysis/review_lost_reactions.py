"""The reactions upstream eQuilibrator can compute and our build cannot.

Two failure modes, and only one of them is a curation signal:

  outside CC span  the reaction still parses and balances, but retraining moved
                   the boundary of the constrained subspace and this reaction
                   fell outside it. Nothing is wrong with the structures.

  unbalanced       the reaction balanced under eQuilibrator's structures and
                   does not under ModelSEED's. That is a direct, self-reporting
                   disagreement about a molecular formula -- exactly what the
                   structure checklist exists to adjudicate.

For the unbalanced ones this prints the atom residual and, per participant,
the formula on each side, so the offending compound is visible without
leaving the terminal.
"""

import os
import csv
import re
import sqlite3
import sys
from pathlib import Path

# ROOT is the eQuilibrator WORKING TREE, not this repository: these analyses read
# caches and fitted parameters that live there and are far too large to commit.
# It was __file__/../.. while this script lived in that tree; now it must be named
# explicitly. Overridable so the analysis is not pinned to one host.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))
sys.path.insert(0, str(ROOT))

OURS = ROOT / "data/modelseed_cache_b/compounds.sqlite"
UP = Path("/home/seaver/.cache/equilibrator/compounds.sqlite")
OUT = ROOT / "data/lost_reactions_review.tsv"


def load(p):
    fh = open(p)
    fh.readline()
    return {r["reaction_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def formula_of(con, acc):
    q = """select c.formula_cache, c.inchi_key from (
             select coalesce(c.atom_bag,'') as formula_cache, c.inchi_key, i.accession
               from compound_identifiers i
               join registries r on r.id = i.registry_id
               join compounds c on c.id = i.compound_id
              where r.namespace like 'seed%') c
            where c.accession = ?"""
    try:
        row = con.execute(q, (acc,)).fetchone()
    except sqlite3.OperationalError:
        return None, None
    return row if row else (None, None)


def main():
    up, new = load(ROOT / "data/modelseed_energies_upstream.tsv"), \
              load(ROOT / "data/modelseed_energies.tsv")
    lost = [k for k in up if up[k]["status"] == "ok"
            and new.get(k, {}).get("status", "").split(":")[0] != "ok"]
    print(f"reactions upstream computes and we do not: {len(lost)}\n")

    from collections import Counter
    by = Counter(new[k]["status"].split(":")[0] for k in lost)
    for k, n in by.most_common():
        print(f"  {k:<22} {n}")

    co, cu = sqlite3.connect(f"file:{OURS}?mode=ro", uri=True), \
             sqlite3.connect(f"file:{UP}?mode=ro", uri=True)

    rows = []
    print("\n" + "=" * 72)
    print("UNBALANCED -- structure disagreements worth curating")
    print("=" * 72)
    for rid in sorted(lost):
        st = new[rid]["status"].split(":")[0]
        if st != "unbalanced":
            rows.append({"reaction_id": rid, "name": new[rid]["name"],
                         "our_status": st, "upstream_dg": up[rid]["dg_prime_kJ_per_mol"],
                         "upstream_sigma": up[rid]["uncertainty_kJ_per_mol"],
                         "differing_compounds": "", "formula": new[rid]["formula"]})
            continue
        accs = re.findall(r"seed:(cpd\d+)", new[rid]["formula"])
        diffs = []
        for a in dict.fromkeys(accs):
            fo, ko = formula_of(co, a)
            fu, ku = formula_of(cu, a)
            if ko and ku and ko != ku:
                diffs.append(f"{a}(ours={ko},up={ku})")
        print(f"\n{rid}  {new[rid]['name'][:58]}")
        print(f"  upstream dG = {up[rid]['dg_prime_kJ_per_mol']} +/- "
              f"{up[rid]['uncertainty_kJ_per_mol']}")
        print(f"  {new[rid]['formula'][:100]}")
        print(f"  differing structures: {', '.join(d.split('(')[0] for d in diffs) or '(none by InChIKey)'}")
        rows.append({"reaction_id": rid, "name": new[rid]["name"],
                     "our_status": st, "upstream_dg": up[rid]["dg_prime_kJ_per_mol"],
                     "upstream_sigma": up[rid]["uncertainty_kJ_per_mol"],
                     "differing_compounds": ";".join(diffs),
                     "formula": new[rid]["formula"]})

    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
