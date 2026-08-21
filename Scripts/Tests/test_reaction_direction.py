#!/usr/bin/env python3
"""Regression test for the Thermodynamics reaction-direction pipeline.

Steps:
1. Bootstrap a read-only baseline by extracting `Biochemistry/reaction_*.json`
   and `compound_*.json` from a reference git branch (default `origin/dev`)
   into `Scripts/Tests/dev_baseline/`. Cached; re-extracted only on demand.
2. Run the Thermodynamics pipeline (same commands as `Rerun_Thermodynamics.sh`).
   This mutates `Biochemistry/*.json` in place — expected behavior.
3. Compare the current tree against the baseline. Checks four things per
   shared reaction:
      * top-level `reversibility`
      * `thermodynamics['Group contribution'][2]` (GC operator)
      * `thermodynamics['eQuilibrator'][2]`       (EQ operator)
      * `thermodynamics['dGPredictor'][2]`        (DGP operator)
   Print a concise summary and exit 0 on PASS, 1 on FAIL.

Usage:
    ./test_reaction_direction.py                    # bootstrap if needed, run, compare
    ./test_reaction_direction.py --no-run           # skip pipeline; just diff
    ./test_reaction_direction.py --refresh-baseline # re-pull baseline
    ./test_reaction_direction.py --baseline-ref origin/master
"""
import argparse
import json
import os
import subprocess
import sys
from glob import glob

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
BIOCHEM_DIR = os.path.join(REPO_ROOT, 'Biochemistry')
THERMO_DIR = os.path.join(REPO_ROOT, 'Scripts', 'Thermodynamics')
BASELINE_DIR = os.path.join(THIS_DIR, 'dev_baseline')

SENTINEL_DG = 10000000

SOURCES = [
    ('GC',  'Group contribution'),
    ('EQ',  'eQuilibrator'),
    ('DGP', 'dGPredictor'),
]

PIPELINE = [
    ['./Update_Compound_GroupContribution_Energies.py'],
    ['./Update_Reaction_GroupContribution_Energies.py'],
    ['./Estimate_Reaction_Reversibility.py', 'GC'],
    ['./Update_Compound_eQuilibrator_Energies.py'],
    ['./Update_Reaction_eQuilibrator_Energies.py'],
    ['./Estimate_Reaction_Reversibility.py', 'EQ'],
    ['./Update_Reaction_dGPredictor_Energies.py'],
    ['./Add_Reaction_Thermodynamics_Operators.py'],
]


def bootstrap_baseline(ref):
    os.makedirs(BASELINE_DIR, exist_ok=True)
    listing = subprocess.run(
        ['git', 'ls-tree', '-r', '--name-only', ref, 'Biochemistry/'],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    names = [n for n in listing.stdout.splitlines()
             if (n.startswith('Biochemistry/reaction_')
                 or n.startswith('Biochemistry/compound_'))
             and n.endswith('.json')]
    if not names:
        sys.exit('ERROR: no reaction_*.json / compound_*.json at ref ' + ref)
    print(f'Bootstrapping baseline from {ref} ({len(names)} files)...')
    for name in names:
        blob = subprocess.run(['git', 'show', f'{ref}:{name}'],
                              cwd=REPO_ROOT, capture_output=True, check=True)
        with open(os.path.join(BASELINE_DIR, os.path.basename(name)), 'wb') as fh:
            fh.write(blob.stdout)


def baseline_present():
    return (bool(glob(os.path.join(BASELINE_DIR, 'reaction_*.json')))
            and bool(glob(os.path.join(BASELINE_DIR, 'compound_*.json'))))


def run_pipeline():
    print(f'Running Thermodynamics pipeline in {THERMO_DIR}...')
    for cmd in PIPELINE:
        print('  $ ' + ' '.join(cmd))
        subprocess.run(cmd, cwd=THERMO_DIR, check=True)


def load_reactions(directory):
    out = {}
    paths = sorted(glob(os.path.join(directory, 'reaction_*.json')))
    if not paths:
        sys.exit('ERROR: no reaction_*.json in ' + directory)
    for path in paths:
        with open(path) as fh:
            for rxn in json.load(fh):
                out[rxn['id']] = rxn
    return out


def source_operator(rxn, label):
    """Return the operator from thermodynamics[label][2], or None if absent
    / sentinel / length-2 legacy entry."""
    thermo = rxn.get('thermodynamics')
    if not isinstance(thermo, dict):
        return None
    sub = thermo.get(label)
    if not sub or sub[0] is None or sub[0] == SENTINEL_DG or len(sub) < 3:
        return None
    return sub[2]


def compare(current, baseline, max_show):
    cur_ids = set(current)
    base_ids = set(baseline)
    only_cur = cur_ids - base_ids
    only_base = base_ids - cur_ids
    shared = cur_ids & base_ids

    rev_mismatch = [rid for rid in sorted(shared)
                    if baseline[rid].get('reversibility')
                    != current[rid].get('reversibility')]

    per_source = {}
    for _, label in SOURCES:
        mismatches = []
        for rid in sorted(shared):
            b_op = source_operator(baseline[rid], label)
            c_op = source_operator(current[rid], label)
            if b_op is not None and c_op is not None and b_op != c_op:
                mismatches.append((rid, b_op, c_op))
        per_source[label] = mismatches

    print('\n' + '=' * 60)
    print('Reaction direction comparison')
    print('=' * 60)
    print(f'  Reactions in baseline               : {len(baseline)}')
    print(f'  Reactions in current                : {len(current)}')
    print(f'  Shared reactions                    : {len(shared)}')
    print(f'  Only in current                     : {len(only_cur)}')
    print(f'  Only in baseline                    : {len(only_base)}')
    print(f'  Top-level reversibility mismatches  : {len(rev_mismatch)}')
    for _, label in SOURCES:
        print(f'  {label:35}: {len(per_source[label]):>4} operator mismatches')

    if rev_mismatch:
        print(f'\nTop-level reversibility mismatches (first {max_show}):')
        print(f'  {"rxn_id":<10} {"baseline":<10} {"current":<10} name')
        for rid in rev_mismatch[:max_show]:
            b_rev = baseline[rid].get('reversibility') or '?'
            c_rev = current[rid].get('reversibility') or '?'
            name = current[rid].get('name', '')[:60]
            print(f'  {rid:<10} {b_rev:<10} {c_rev:<10} {name}')

    for _, label in SOURCES:
        mm = per_source[label]
        if mm:
            print(f'\n{label} operator mismatches (first {max_show}):')
            for rid, b_op, c_op in mm[:max_show]:
                print(f'  {rid:<10} {b_op} → {c_op}')

    if only_cur:
        print(f'\nOnly in current (first {max_show}):')
        for rid in sorted(only_cur)[:max_show]:
            print(f'  + {rid}')
    if only_base:
        print(f'\nOnly in baseline (first {max_show}):')
        for rid in sorted(only_base)[:max_show]:
            print(f'  - {rid}')

    ok = (not rev_mismatch and not only_cur and not only_base
          and all(not per_source[label] for _, label in SOURCES))
    print('\nRESULT: ' + ('PASS' if ok else 'FAIL'))
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--baseline-ref', default='origin/dev',
                   help='git ref to pull baseline from (default: origin/dev)')
    p.add_argument('--refresh-baseline', action='store_true',
                   help='re-extract the baseline')
    p.add_argument('--no-run', action='store_true',
                   help='skip the pipeline; just diff')
    p.add_argument('--display-num', type=int, default=20, metavar='N',
                   help='max rows per diff section (default: 20)')
    args = p.parse_args()

    if args.display_num < 0:
        p.error('--display-num must be non-negative')

    if args.refresh_baseline or not baseline_present():
        bootstrap_baseline(args.baseline_ref)
    else:
        print(f'Baseline present at {BASELINE_DIR} '
              f'(use --refresh-baseline to rebuild)')

    if not args.no_run:
        run_pipeline()
    else:
        print('Skipping pipeline run (--no-run)')

    baseline = load_reactions(BASELINE_DIR)
    current = load_reactions(BIOCHEM_DIR)
    ok = compare(current, baseline, max_show=args.display_num)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
