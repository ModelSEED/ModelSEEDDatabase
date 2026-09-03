#!/usr/bin/env python
"""Rewrite each compound's per-source ``pkas`` dict into the index-free form.

Two changes, both driven by measurement rather than taste.

**Atom indices are dropped.** Nothing in the energy path consumes one:
eQuilibrator stores dissociation constants as a positional list and only takes
its length, compares entries to the pH and sums a slice; equilibrator-assets
takes SMILES -> list of floats; dGPredictor does not use pKas; Group
Contribution does not reference them. The one reader was
``Build_Biochemistry_Object.py``, and 64.2% of the atom keys it published
pointed at carbon, phosphorus, or an index past the end of the molecule --
because Marvin reorders atoms on import, so its indices never described our
structures. The indices are not destroyed, only left where they mean
something: the per-source files under ``Biochemistry/Structures/<source>/pkas/``,
each in the atom space it was computed in.

**Each source declares a ``kind``.** ``microscopic`` for per-site predictions
on one protonation state, which is what Marvin and MolGpKa produce;
``macroscopic`` for sequential dissociation constants of the whole molecule,
which is what measured values are and what eQuilibrator's transform actually
requires. eQuilibrator reads any list as a macroscopic ladder, so handing it
microscopic values is a category error rather than an inaccuracy -- it moved
one LPS reaction by 355 kcal/mol.

The served top-level ``pka`` / ``pkb`` strings are NOT touched, so the KBase
biochemistry object is unchanged by this script. Correcting that is a separate
decision with downstream consumers (see the audit in the paper notes).

Idempotent: a record already in two-field form is left alone.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pka_encoding import MACROSCOPIC, MICROSCOPIC, migrate

# Every source we currently hold predicts per-site values on one protonation
# state. A macroscopic source (measured literature values) will set this
# explicitly when it is added.
KINDS = {"Marvin": MICROSCOPIC, "MolGpKa": MICROSCOPIC, "Literature": MACROSCOPIC}

BIOCHEM = Path(__file__).resolve().parents[2] / "Biochemistry"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    seen = {"records": 0, "already": 0, "migrated": 0, "unknown_source": 0}
    for path in sorted(BIOCHEM.glob("compound_*.json")):
        entries = json.load(path.open())
        touched = False
        for e in entries:
            pkas = e.get("pkas")
            if not isinstance(pkas, dict):
                continue
            for source, rec in pkas.items():
                if not isinstance(rec, dict):
                    continue
                seen["records"] += 1
                kind = KINDS.get(source)
                if kind is None:
                    seen["unknown_source"] += 1
                    continue
                before = (rec.get("pKa"), rec.get("pKb"), rec.get("kind"))
                rec["kind"] = kind
                for field in ("pKa", "pKb"):
                    rec[field] = migrate(rec.get(field) or "")
                after = (rec.get("pKa"), rec.get("pKb"), rec.get("kind"))
                if before == after:
                    seen["already"] += 1
                else:
                    seen["migrated"] += 1
                    touched = True
        if touched and not a.dry_run:
            # byte-for-byte the format BiochemPy.Compounds.saveCompounds
            # writes: indent 4, sorted keys, and NO trailing newline. Anything
            # else reformats all 46 files and buries the real change.
            with path.open("w") as fh:
                fh.write(json.dumps(entries, indent=4, sort_keys=True))

    for k, v in seen.items():
        print(f"  {k:<16}{v:>9,}")
    if a.dry_run:
        print("\n(dry run - nothing written)")


if __name__ == "__main__":
    main()
