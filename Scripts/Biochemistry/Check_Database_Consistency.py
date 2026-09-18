#!/usr/bin/env python3
"""Verify that derived fields in the biochemistry database are consistent
with their source-of-truth fields. Read-only -- does not modify any files.

Source-of-truth fields are values you change deliberately (reaction
`reversibility`; SMILES/InChI/InChIKey in
`Biochemistry/Structures/Unique_ModelSEED_Structures.txt`; etc.).

Derived fields are STORED copies that should match what the canonical
regeneration code (Reprint_Biochemistry.py +
BiochemPy.Reactions.rebuildReaction, BiochemPy.Compounds.saveCompounds)
would produce. They drift when a source field is updated without
running the corresponding reprint/regeneration step.

Invariants checked:

  [reactions]
    The `equation` and `definition` fields' direction arrow must match
    what rebuildReaction derives from `reversibility`:
        reversibility='<'  -> '<='
        reversibility='>'  -> '=>'
        anything else      -> '<=>'
    Fix: run Scripts/Biochemistry/Reprint_Biochemistry.py

  [compounds]
    The `smiles` and `inchikey` fields must match the values in
    Biochemistry/Structures/Unique_ModelSEED_Structures.txt (as loaded
    by BiochemPy.Compounds.loadStructures for the ModelSEED source).
    Fix: run Scripts/Structures/Update_Compound_Structures_Formulas_Charge.py

Exit codes:
  0   all invariants pass
  1   one or more invariants violated (drift detected)
  2   internal error (file missing, parse failure, etc.)

Usage:
  python3 Scripts/Biochemistry/Check_Database_Consistency.py
  python3 Scripts/Biochemistry/Check_Database_Consistency.py --only reactions
  python3 Scripts/Biochemistry/Check_Database_Consistency.py --only compounds
  python3 Scripts/Biochemistry/Check_Database_Consistency.py --max-report 100

Designed for use locally before pushing, in a pre-push git hook, and
in CI. Runs in a few seconds against the full database.
"""

import argparse
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(
    os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python')))

try:
    from BiochemPy import Compounds, Reactions
except ImportError as e:
    print(f"ERROR: cannot import BiochemPy ({e}). "
          f"Ensure Libs/Python is on PYTHONPATH or run this script from "
          f"its installed location in Scripts/Biochemistry/.", file=sys.stderr)
    sys.exit(2)


# ────────────────────────────────────────────────────────────────────
# Reaction invariant
# ────────────────────────────────────────────────────────────────────

def derive_arrow(reversibility):
    """Mirror BiochemPy.Reactions.rebuildReaction direction logic."""
    if reversibility == '<':
        return '<='
    if reversibility == '>':
        return '=>'
    return '<=>'


def stored_arrow(equation):
    """Extract the direction arrow present in an equation/definition string.
    Order matters: must check '<=>' before its substrings '=>' and '<='."""
    if not equation:
        return '?'
    if '<=>' in equation:
        return '<=>'
    if '=>' in equation:
        return '=>'
    if '<=' in equation:
        return '<='
    return '?'


def check_reactions(max_report):
    """Verify reaction equation/definition arrows match reversibility.
    Returns the number of reactions found to be in drift."""
    print("[reactions] loading...")
    try:
        rxns = Reactions().loadReactions()
    except Exception as e:
        print(f"  ERROR loading reactions: {e}", file=sys.stderr)
        return -1

    drift = []
    for rid, r in rxns.items():
        rev = r.get('reversibility')
        expected = derive_arrow(rev)
        # Check both rendered fields; report a reaction once even if both drift.
        for field in ('equation', 'definition'):
            actual = stored_arrow(r.get(field, ''))
            if actual != expected:
                drift.append((rid, field, rev, actual, expected))
                break

    if not drift:
        print(f"[reactions] OK ({len(rxns)} reactions checked)")
        return 0

    print(f"[reactions] DRIFT: {len(drift)} of {len(rxns)} reactions have "
          f"stale equation/definition arrows")
    print(f"  Fix: python3 Scripts/Biochemistry/Reprint_Biochemistry.py")
    print(f"  Examples (up to {max_report}):")
    for rid, field, rev, actual, expected in drift[:max_report]:
        print(f"    {rid}: reversibility={rev!r}  "
              f"stored {field} arrow={actual}  expected={expected}")
    if len(drift) > max_report:
        print(f"    ... and {len(drift) - max_report} more")
    return len(drift)


# ────────────────────────────────────────────────────────────────────
# Compound invariant
# ────────────────────────────────────────────────────────────────────

def check_compounds(max_report):
    """Verify compound smiles/inchikey fields match Unique_ModelSEED_Structures.txt.
    Returns the number of compounds found to be in drift."""
    print("[compounds] loading...")
    try:
        helper = Compounds()
        compounds = helper.loadCompounds()
        structures = helper.loadStructures(["SMILE", "InChIKey"], ["ModelSEED"])
    except Exception as e:
        print(f"  ERROR loading compounds/structures: {e}", file=sys.stderr)
        return -1

    drift = []
    for cpd, c in compounds.items():
        expected_smiles = ''
        expected_key = ''
        if cpd in structures:
            if 'SMILE' in structures[cpd]:
                expected_smiles = sorted(structures[cpd]['SMILE'].keys())[0]
            if 'InChIKey' in structures[cpd]:
                expected_key = sorted(structures[cpd]['InChIKey'].keys())[0]

        stored_smiles = c.get('smiles', '') or ''
        stored_key = c.get('inchikey', '') or ''

        if stored_smiles != expected_smiles or stored_key != expected_key:
            drift.append((cpd, stored_smiles, expected_smiles,
                          stored_key, expected_key))

    if not drift:
        print(f"[compounds] OK ({len(compounds)} compounds checked)")
        return 0

    print(f"[compounds] DRIFT: {len(drift)} of {len(compounds)} compounds "
          f"have stale smiles or inchikey")
    print(f"  Fix: python3 Scripts/Structures/"
          f"Update_Compound_Structures_Formulas_Charge.py")
    print(f"  Examples (up to {max_report}):")
    for cpd, ss, es, sk, ek in drift[:max_report]:
        diffs = []
        if ss != es:
            s_display = ss[:50] + ('...' if len(ss) > 50 else '')
            e_display = es[:50] + ('...' if len(es) > 50 else '')
            diffs.append(f"smiles: stored={s_display!r} expected={e_display!r}")
        if sk != ek:
            diffs.append(f"inchikey: stored={sk!r} expected={ek!r}")
        print(f"    {cpd}: {'; '.join(diffs)}")
    if len(drift) > max_report:
        print(f"    ... and {len(drift) - max_report} more")
    return len(drift)


# ────────────────────────────────────────────────────────────────────
# Entry point
# ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Check biochemistry database for derived-field drift",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--only', choices=['reactions', 'compounds', 'all'], default='all',
        help='Restrict to one invariant family (default: all)')
    parser.add_argument(
        '--max-report', type=int, default=20,
        help='Maximum number of example violations to print per family '
             '(default: 20)')
    args = parser.parse_args()

    print(f"Checking database at {os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))}\n")

    total_drift = 0
    if args.only in ('reactions', 'all'):
        n = check_reactions(args.max_report)
        if n < 0:
            return 2
        total_drift += n
        print()
    if args.only in ('compounds', 'all'):
        n = check_compounds(args.max_report)
        if n < 0:
            return 2
        total_drift += n
        print()

    if total_drift == 0:
        print("All consistency invariants pass.")
        return 0

    print(f"FAIL: {total_drift} total drift instance(s) detected.")
    return 1


if __name__ == '__main__':
    sys.exit(main())
