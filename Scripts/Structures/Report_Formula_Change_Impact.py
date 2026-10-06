#!/usr/bin/env python3
"""Report the mass-balance impact of proposed compound formula/charge changes
on every reaction that uses any of the affected compounds.

Given a TSV of proposed changes (schema compatible with
acps_formula_charge.tsv: at minimum an `ID` (or `cpd_id`) column,
`formula`, and `charge`), the script:

  1. Loads current compound formulas from Biochemistry/compound_*.json
     via BiochemPy.Compounds.
  2. Builds an "after" formula/charge lookup by overlaying the proposed
     changes.
  3. For every reaction that involves any proposed compound, computes
     mass and charge balance BEFORE (current formulas) and AFTER
     (proposed formulas).
  4. Buckets each reaction as one of:

       still_balanced      balanced -> balanced             (no impact)
       newly_balanced      imbalanced -> balanced           (fix!)
       newly_imbalanced    balanced -> imbalanced           (regression!)
       still_imbalanced    imbalanced -> imbalanced         (no help)

  5. Emits a report to stdout (or --output) with per-bucket details and
     a summary of counts.

Wildcards ('R' in formulas) are treated as an element and must balance
the same as any other atom. Reactions where both sides have equal R
counts stay balanced; reactions where an ACP override changes R count
will surface as newly-imbalanced.

Usage:
  python3 Scripts/Structures/Report_Formula_Change_Impact.py proposed.tsv
  python3 Scripts/Structures/Report_Formula_Change_Impact.py \\
      Biochemistry/Curation/overrides/acps_formula_charge.tsv \\
      --output impact.txt
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(
    os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python')))

try:
    from BiochemPy import Compounds, Reactions
except ImportError as e:
    print(f"ERROR: cannot import BiochemPy ({e}). Run from a location "
          f"where Libs/Python/BiochemPy is importable.", file=sys.stderr)
    sys.exit(2)


def parse_formula(f):
    """Return {element: count} from a formula string like 'C11H21N2O7PRS'.
    Treats 'R' as a first-class element (wildcard placeholder)."""
    if not f or f in ('null', 'NULL'):
        return {}
    out = {}
    for m in re.finditer(r'([A-Z][a-z]?|R|\*)(\d*)', f):
        elem = m.group(1)
        n = int(m.group(2)) if m.group(2) else 1
        out[elem] = out.get(elem, 0) + n
    return out


def balance_delta(reagents, formula_lookup, charge_lookup):
    """Compute (atom_delta_dict, charge_delta) for a reaction.
    Balanced means all atom deltas and charge delta are zero."""
    atoms = defaultdict(float)
    charge = 0.0
    for rgt in reagents:
        cpd = rgt.get('compound', '')
        try:
            coef = float(rgt.get('coefficient', 0))
        except (TypeError, ValueError):
            coef = 0.0
        f = formula_lookup.get(cpd, '')
        c = charge_lookup.get(cpd, 0)
        for elem, n in parse_formula(f).items():
            atoms[elem] += coef * n
        try:
            charge += coef * float(c)
        except (TypeError, ValueError):
            pass
    return {k: v for k, v in atoms.items() if abs(v) > 1e-9}, charge


def format_delta(atoms, charge):
    """Compact human-readable delta string; 'OK' if balanced."""
    if not atoms and abs(charge) < 1e-9:
        return 'OK'
    parts = []
    if atoms:
        parts.append(','.join(
            f"{k}:{int(v) if v == int(v) else v:+g}"
            for k, v in sorted(atoms.items())))
    if abs(charge) > 1e-9:
        parts.append(f"chg:{int(charge) if charge == int(charge) else charge:+g}")
    return ' | '.join(parts)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split('\n\n')[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument('proposed',
                        help='TSV of proposed changes '
                             '(cols: ID or cpd_id, formula, charge; header required)')
    parser.add_argument('--output', default='-',
                        help='Report output path (default: stdout)')
    parser.add_argument('--max-per-bucket', type=int, default=50,
                        help='Max reactions to detail per bucket (default: 50)')
    args = parser.parse_args()

    # Load proposed changes
    proposed = {}
    with open(args.proposed) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        id_col = 'ID' if 'ID' in reader.fieldnames else 'cpd_id'
        for row in reader:
            cpd = row.get(id_col, '').strip()
            if not cpd:
                continue
            f = row.get('formula', '').strip()
            c = row.get('charge', '').strip()
            if f and c:
                proposed[cpd] = (f, c)
    print(f"Loaded {len(proposed)} proposed compound formula changes",
          file=sys.stderr)

    if not proposed:
        print("No valid proposed changes; exiting.", file=sys.stderr)
        return 1

    # Load current compound formulas
    print("Loading current compounds...", file=sys.stderr)
    ch = Compounds()
    compounds = ch.loadCompounds()
    formula_before = {c: (compounds[c].get('formula') or '') for c in compounds}
    charge_before = {c: (compounds[c].get('charge') or 0) for c in compounds}

    # Build 'after' by overlaying proposed changes
    formula_after = dict(formula_before)
    charge_after = dict(charge_before)
    for cpd, (f, c) in proposed.items():
        formula_after[cpd] = f
        try:
            charge_after[cpd] = int(c)
        except ValueError:
            charge_after[cpd] = c

    # Load reactions and find affected ones
    print("Loading reactions...", file=sys.stderr)
    rh = Reactions()
    reactions = rh.loadReactions()
    affected = []
    proposed_set = set(proposed)
    for rxn_id, r in reactions.items():
        cpd_ids = set(str(r.get('compound_ids', '')).split(';'))
        if cpd_ids & proposed_set:
            affected.append((rxn_id, r))
    print(f"Reactions involving any proposed compound: {len(affected)} "
          f"of {len(reactions)}", file=sys.stderr)

    # Bucket by before/after balance
    buckets = defaultdict(list)
    for rxn_id, r in affected:
        rgts = r.get('stoichiometry', []) or []
        atoms_b, chg_b = balance_delta(rgts, formula_before, charge_before)
        atoms_a, chg_a = balance_delta(rgts, formula_after, charge_after)
        balanced_b = not atoms_b and abs(chg_b) < 1e-9
        balanced_a = not atoms_a and abs(chg_a) < 1e-9
        involved = sorted(
            proposed_set & {rg.get('compound', '') for rg in rgts})
        entry = {
            'rxn_id': rxn_id,
            'involved': ';'.join(involved),
            'before': format_delta(atoms_b, chg_b),
            'after': format_delta(atoms_a, chg_a),
        }
        if balanced_b and balanced_a:
            b = 'still_balanced'
        elif not balanced_b and balanced_a:
            b = 'newly_balanced'
        elif balanced_b and not balanced_a:
            b = 'newly_imbalanced'
        else:
            b = 'still_imbalanced'
        buckets[b].append(entry)

    # Emit report
    out = sys.stdout if args.output == '-' else open(args.output, 'w')
    for bucket in ['newly_imbalanced', 'newly_balanced',
                   'still_imbalanced', 'still_balanced']:
        entries = buckets.get(bucket, [])
        print(f"\n=== {bucket} ({len(entries)}) ===", file=out)
        for e in entries[:args.max_per_bucket]:
            print(f"  {e['rxn_id']:<10} cpds={e['involved']:<25} "
                  f"BEFORE: {e['before']:<30} AFTER: {e['after']}", file=out)
        if len(entries) > args.max_per_bucket:
            print(f"  ... and {len(entries) - args.max_per_bucket} more",
                  file=out)
    if out is not sys.stdout:
        out.close()

    print(f"\nSummary:", file=sys.stderr)
    for k in ('newly_imbalanced', 'newly_balanced',
              'still_imbalanced', 'still_balanced'):
        print(f"  {k}: {len(buckets.get(k, []))}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
