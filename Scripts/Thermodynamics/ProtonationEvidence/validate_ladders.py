#!/usr/bin/env python
"""Check ``literature_ladders.tsv`` against the magnesium guard, and refuse
uncited rows.

Two independent checks, because the table can fail in two unrelated ways.

1. **Citation.** Every row must carry a non-empty ``citation``. This is the
   whole point of the file: a table where some rows are cited and others are
   remembered is worse than no table, because a reviewer cannot tell which is
   which. An uncited row is an error, not a warning.

2. **Guard admissibility.** A ladder can be perfectly well cited and still be
   useless, because ``build_modelseed_cache.py`` will refuse to install it.
   The guard is reproduced here exactly::

       major     = sum(1 for v in ladder if v > 7.0)
       reachable = {i - major + n_h for i in range(len(ladder) + 1)}
       proceed if set(mg_protons) <= reachable

   Equivalently the ladder needs ``max(0, n_h - min(mg))`` values above pH 7
   and ``max(0, max(mg) - n_h)`` at or below it. A compound whose ladder fails
   keeps its ChemAxon protonation, so FAIL here means the curation effort for
   that compound bought nothing yet.

Ladders are assembled per ``(seed_id, source)`` and never merged across
sources: mixing a value measured at I = 0 with one measured in 0.1 M KCl
produces a ladder that is neither, and the resolver takes one source per
compound anyway.

Exit status is non-zero if any row is uncited or any assembled ladder that has
magnesium data fails the guard, so this is usable as a pre-commit gate.
"""
import csv
import pickle
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

from paths import EQ, THERMO_CACHE, DATA_OUT, require

CACHE = THERMO_CACHE
TABLE = DATA_OUT / "literature_ladders.tsv"

sys.path.insert(0, str(EQ / "tools"))
from modelseed_pkas import cache_seed_identifiers  # noqa: E402


def compound_facts(seed_ids):
    """seed_id -> (n_h, sorted mg proton counts) from the cache."""
    sm = cache_seed_identifiers(str(CACHE))
    con = sqlite3.connect(str(CACHE))
    facts = {}
    for sid in seed_ids:
        ent = sm.get(sid)
        if not ent:
            continue
        cid = ent[0]
        row = con.execute(
            "select atom_bag from compounds where id=?", (cid,)).fetchone()
        if not row:
            continue
        bag = row[0]
        bag = pickle.loads(bag) if isinstance(bag, (bytes, bytearray)) else (bag or {})
        n_h = bag.get("H", 0) if isinstance(bag, dict) else 0
        mg = sorted({x[0] for x in con.execute(
            "select number_protons from magnesium_dissociation_constant "
            "where compound_id=?", (cid,))})
        facts[sid] = (n_h, mg)
    return facts


def guard(ladder, n_h, mg):
    major = sum(1 for v in ladder if v > 7.0)
    reachable = {i - major + n_h for i in range(len(ladder) + 1)}
    return set(mg) <= reachable


def main():
    if not TABLE.exists():
        sys.exit(f"missing {TABLE}")
    rows = list(csv.DictReader(TABLE.open(), delimiter="\t"))
    uncited = [r for r in rows if not (r.get("citation") or "").strip()]

    ladders = defaultdict(list)
    names = {}
    for r in rows:
        try:
            ladders[(r["seed_id"], r["source"])].append(float(r["pka_value"]))
        except (TypeError, ValueError):
            uncited.append(r)          # unparseable value is equally unusable
            continue
        names[r["seed_id"]] = r["name"]

    facts = compound_facts({s for s, _ in ladders})

    print(f"rows {len(rows)}   compounds {len(names)}   "
          f"(compound, source) ladders {len(ladders)}")
    print(f"uncited or unparseable rows: {len(uncited)}")
    for r in uncited[:10]:
        print(f"  UNCITED  {r.get('seed_id')} {r.get('step')} {r.get('pka_value')}")
    print()

    failures = 0
    hdr = f"{'compound':10s} {'source':9s} {'name':26s} {'ladder':32s} {'nH':>3} {'mg':16s} verdict"
    print(hdr)
    print("-" * len(hdr))
    for (sid, src), vals in sorted(ladders.items()):
        vals = sorted(vals, reverse=True)
        n_h, mg = facts.get(sid, (None, None))
        if n_h is None:
            verdict = "not in cache"
        elif not mg:
            verdict = "n/a (no Mg data)"
        elif guard(vals, n_h, mg):
            verdict = "PASS"
        else:
            need_above = max(0, n_h - min(mg))
            need_below = max(0, max(mg) - n_h)
            have_above = sum(1 for v in vals if v > 7.0)
            verdict = (f"FAIL need >={need_above} above / >={need_below} below pH7, "
                       f"have {have_above} / {len(vals) - have_above}")
            failures += 1
        print(f"{sid:10s} {src:9s} {names[sid][:26]:26s} "
              f"{' '.join(f'{v:.2f}' for v in vals)[:32]:32s} "
              f"{n_h if n_h is not None else '-':>3} "
              f"{(';'.join(map(str, mg)) if mg else '-')[:16]:16s} {verdict}")

    print()
    print(f"guard failures: {failures}")
    if uncited or failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
