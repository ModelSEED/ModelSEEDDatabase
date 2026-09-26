#!/usr/bin/env python
"""Replace the MolGpKa entry in each compound's ``pkas`` dict from a run
registered in ``Biochemistry/Structures/sources.yaml``.

Replace, not merge. The stored entry came from an earlier run with different
coverage (22,395 compounds against 32,201) and a different site convention --
it recorded the ionizable HYDROGEN where the current run records the heavy atom
bearing it. Merging two runs that disagree about what an index means would
produce a record that is true of neither.

Values are written in the index-free two-field encoding, because nothing in the
energy path consumes an atom index; the indices stay in the source file under
``Biochemistry/Structures/ModelSEED/pkas/``, which is the only place they mean
anything. See ``pka_encoding`` for why.

Sites are kept over -2..16 rather than 0..14. eQuilibrator counts sites above
the reported pH to place the major microspecies, so a group at pKa 14.9 is
protonated at pH 7 and carries a proton; discarding it moves dG'0 by a full
RT ln(10) x pH = 9.55 kcal/mol. Multiple protons on one heavy atom collapse to
a single entry -- an NH2 has two hydrogens and one first dissociation -- which
is what the heavy-atom remap in the source run already encodes.
"""
import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pka_encoding import MICROSCOPIC, encode

ROOT = Path(__file__).resolve().parents[2]
BIOCHEM = ROOT / "Biochemistry"
STRUCT = BIOCHEM / "Structures"
LABEL = "MolGpKa"
LO, HI = -2.0, 16.0


def load_run():
    """The ModelSEED pKa file flagged consumed_by_production in sources.yaml."""
    cfg = yaml.safe_load((STRUCT / "sources.yaml").open())
    entries = [e for e in cfg["sources"]["ModelSEED"]["pkas"]
               if e.get("consumed_by_production")]
    if len(entries) != 1:
        raise SystemExit(
            f"expected exactly one production ModelSEED pKa file, found "
            f"{len(entries)}: {[e['file'] for e in entries]}")
    entry = entries[0]
    path = STRUCT / "ModelSEED" / entry["file"]
    print(f"source: {entry['file']}  ({entry['tool']}/{entry['version']})")

    out = defaultdict(lambda: {"pKa": [], "pKb": []})
    dropped = kept = 0
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            for tok in row["value"].split(";"):
                if not tok:
                    continue
                value = float(tok.rsplit(":", 1)[-1])
                if LO <= value <= HI:
                    out[row["external_id"]][row["kind"]].append(value)
                    kept += 1
                else:
                    dropped += 1
    print(f"  sites kept {kept:,}   outside {LO}..{HI} {dropped:,}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    run = load_run()
    seen = {"compounds with a run entry": 0, "replaced": 0,
            "added": 0, "removed (not in run)": 0, "unchanged": 0}

    for path in sorted(BIOCHEM.glob("compound_*.json")):
        entries = json.load(path.open())
        touched = False
        for e in entries:
            pkas = e.setdefault("pkas", {}) if e["id"] in run else e.get("pkas")
            if not isinstance(pkas, dict):
                continue
            if e["id"] in run:
                seen["compounds with a run entry"] += 1
                rec = {"kind": MICROSCOPIC}
                for field in ("pKa", "pKb"):
                    vals = sorted(run[e["id"]][field], reverse=True)
                    rec[field] = encode(vals) if vals else ""
                before = pkas.get(LABEL)
                if before == rec:
                    seen["unchanged"] += 1
                    continue
                seen["replaced" if before else "added"] += 1
                pkas[LABEL] = rec
                touched = True
            elif LABEL in pkas:
                # the run covered every structure, so absence is a positive
                # statement: this compound has no structure to predict from
                del pkas[LABEL]
                seen["removed (not in run)"] += 1
                touched = True
            if not pkas:
                e.pop("pkas", None)
        if touched and not a.dry_run:
            with path.open("w") as fh:
                fh.write(json.dumps(entries, indent=4, sort_keys=True))

    for k, v in seen.items():
        print(f"  {k:<28}{v:>9,}")
    if a.dry_run:
        print("\n(dry run - nothing written)")


if __name__ == "__main__":
    main()
