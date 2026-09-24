#!/usr/bin/env python
"""Replace protonation rows that break chemistry with their unprotonated source.

THE RULE. A protonation moves protons: between a source structure and its
pH-adjusted form, dH must equal dcharge and no other element may change
(Validate_Protonations.py). A row that fails this did not undergo a
protonation -- in every case found so far Marvin was handed a molecule that
InChI had already taken apart, or a bare S/Se/P/O atom, and protonated the
pieces as free ions. The protonation state of a shattered fragment set is not
defined, so there is no correct protonated row to write. What IS defined is
the source: a self-consistent description at whatever protonation state the
source database chose. So:

    a row that fails the invariant is replaced by its own representation's
    source row -- structure, formula and charge -- and, when that row is the
    InChI, the compound's InChIKey row is re-hashed from it. Nothing else
    changes.

WHY THIS IS SAFE. The rule only touches rows the invariant rejects. Rows where
Marvin legitimately moved protons pass it (dH == dcharge) and are left alone;
that includes the 297 disconnected-InChI compounds in the 26.1 bundle whose
InChI row picks up a /p-2 layer -- the chlorophylls and cobalamins that the
run script's "protonate the InChI form rather than pass it through" decision
was made to protect. Both concerns are met at once: protonate everything, and
keep the result only where it is a protonation.

WHAT IT DOES NOT FIX. A row that passes but reproduces a source that is itself
inconsistent -- MetaCyc CPD-7's SMILES carries [SH] on bridging sulfides that
its InChI does not -- is untouched, because the source is the reference. The
validator's `source` check lists those; they are curation, not repair.

The same function is the gate the run script should apply as it writes each
row, so a future regeneration never ships such a row in the first place;
applying it after the fact is deterministic and needs no Marvin licence.

USAGE
    Repair_Protonation_Rows.py                        # dry run: report only
    Repair_Protonation_Rows.py --out /tmp/repaired    # write repaired tree there
    Repair_Protonation_Rows.py --write                # repair in place; every row
                                                      # replaced is recorded in
                                                      # _reports/<bundle>_passthrough_<source>.tsv
"""
import argparse
import csv
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Validate_Protonations as V   # noqa: E402

BUNDLE_COLS = ["external_id", "type", "structure", "formula", "charge",
               "tool", "tool_version", "ph", "generated_on"]


def repair_bundle(src, bundle, root):
    """Return (rows_out, replaced) for one source's bundle under `root`."""
    V.STRUCT = root
    base = V.load_source(src)
    path = os.path.join(root, src, "protonations", bundle + ".tsv")
    if not os.path.isfile(path):
        return None, []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        rows = list(rd)
    replaced = []
    rekey = {}          # external_id -> key of the source InChI now written
    for r in rows:
        typ = r["type"]
        if typ not in V.ORIGINAL or not r["formula"]:
            continue
        k = (r["external_id"], typ)
        if k not in base:
            continue
        f0, q0, s0 = base[k]
        kind, dH, dq = V.compare(f0, q0, r["formula"], r["charge"])
        if kind == "ok" or kind in V.INFO_KINDS:
            continue
        replaced.append(dict(external_id=r["external_id"], type=typ, kind=kind, dH=dH, dq=dq,
                             was_structure=r["structure"], was_formula=r["formula"], was_charge=r["charge"],
                             now_structure=s0, now_formula=f0, now_charge=q0))
        r["structure"], r["formula"], r["charge"] = s0, f0, q0
        if typ == "InChI":
            rekey[r["external_id"]] = V.inchikey_of(s0)
    # The InChIKey row is a hash of the InChI row and must follow it. The
    # first version of this script did not do this, and elemental sulfur
    # shipped keyed as hydrosulfide -- the key of the protonated InChI that
    # had just been discarded. The gate in Run_Marvin_Protonations.py
    # re-hashes for the same reason.
    for r in rows:
        if r["type"] == "InChIKey" and r["external_id"] in rekey and rekey[r["external_id"]]:
            key = rekey[r["external_id"]]
            if r["structure"] != key:
                replaced.append(dict(external_id=r["external_id"], type="InChIKey",
                                     kind="rehashed from the passed-through InChI", dH="", dq="",
                                     was_structure=r["structure"], was_formula="", was_charge="",
                                     now_structure=key, now_formula="", now_charge=""))
                r["structure"] = key
    return (cols, rows), replaced


def write_bundle(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", default="marvin_26.1_ph7")
    ap.add_argument("--sources", nargs="*", default=list(V.SOURCES))
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--out", help="write a repaired copy of Biochemistry/Structures here")
    g.add_argument("--write", action="store_true", help="repair the bundle in place")
    a = ap.parse_args()

    src_root = V.STRUCT
    if a.out:
        for s in a.sources:
            d = os.path.join(src_root, s)
            if os.path.isdir(d):
                shutil.copytree(d, os.path.join(a.out, s), dirs_exist_ok=True)
        root = a.out
    else:
        root = src_root

    total = 0
    print(f"bundle {a.bundle}   mode: {'write in place' if a.write else 'write to ' + a.out if a.out else 'DRY RUN'}")
    print(f"{'source':<9}{'rows':>8}{'replaced':>10}   by kind")
    for s in a.sources:
        out, replaced = repair_bundle(s, a.bundle, root)
        if out is None:
            continue
        cols, rows = out
        kinds = {}
        for r in replaced:
            kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
        print(f"{s:<9}{len(rows):>8,}{len(replaced):>10}   " + ", ".join(f"{k}={v}" for k, v in sorted(kinds.items())))
        total += len(replaced)
        if (a.out or a.write) and replaced:
            path = os.path.join(root, s, "protonations", a.bundle + ".tsv")
            write_bundle(path, cols, rows)
            # NOT under protonations/: BiochemPy.loadStructures globs
            # protonations/*.tsv as bundles, and a sidecar there would be read
            # as one. Provenance goes with the other run reports instead.
            rep = os.path.join(root, "_reports")
            os.makedirs(rep, exist_ok=True)
            side = os.path.join(rep, f"{a.bundle}_passthrough_{s}.tsv")
            with open(side, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(replaced[0].keys()), delimiter="\t", lineterminator="\n")
                w.writeheader(); w.writerows(replaced)
    print(f"\n{total} rows {'replaced' if (a.out or a.write) else 'would be replaced'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
