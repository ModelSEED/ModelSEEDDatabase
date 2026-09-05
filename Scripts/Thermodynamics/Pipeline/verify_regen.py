#!/usr/bin/env python3
"""Compare a regenerated eQuilibrator table against the one committed in MSD.

Answers the three questions that decide whether the regeneration is safe to
ingest:

  1. Did anything REGRESS -- was 'ok', now is not? Any hit is a bug, not a
     policy change, because the regeneration only ever splits a decline into
     two finer declines. Exception: rows whose old value came from the
     residual branch, which never had a real estimate to lose.
  2. Of the rows that newly carry a value, how many are genuinely usable
     versus effectively unconstrained (sigma at or beyond RMSE_inf)?
  3. Do the rows present in both agree?
"""
import sys
import os
import argparse, csv, sys
from collections import Counter
from pathlib import Path

# The eQuilibrator working tree: caches and fitted parameters, gigabytes, not
# in this repository. Named by environment variable since relocation.
ROOT = Path(os.environ.get("EQUILIBRATOR_DIR",
                           "/scratch/seaver/Claude_Projects/eQuilibrator"))


KCAL = 4.184
RMSE_INF_KCAL = 100000.0 / KCAL     # 23,901

def load(path, idcol):
    out = {}
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith('#'):
            fh.seek(0)
        for r in csv.DictReader(fh, delimiter='\t'):
            out[r[idcol]] = r
    return out

def num(r, col):
    v = (r.get(col) or '').strip()
    if v in ('', 'inf', 'nan'):
        return None
    try:
        return float(v)
    except ValueError:
        return None

def median(xs):
    xs = sorted(xs)
    return xs[len(xs)//2] if xs else float('nan')

def report(old, new, idcol, dgcol, errcol, label):
    print('=' * 72)
    print(label)
    print('=' * 72)
    common = [k for k in new if k in old]
    o_ok = {k for k in common if old[k]['status'] == 'ok'}
    n_ok = {k for k in common if new[k]['status'] == 'ok'}
    print(f'rows compared        : {len(common)}')
    print(f'  ok before          : {len(o_ok)}')
    print(f'  ok after           : {len(n_ok)}   ({len(n_ok)-len(o_ok):+d})')

    # 1. regressions
    reg = sorted(o_ok - n_ok)
    print(f'\nREGRESSIONS (was ok, now not): {len(reg)}')
    for k in reg[:10]:
        print(f'    {k}  -> {new[k]["status"][:58]}')
    if len(reg) > 10:
        print(f'    ... and {len(reg)-10} more')

    # 2. newly valued, split by usability
    gained = sorted(n_ok - o_ok)
    # Bands, not a single arbitrary cut. An earlier version of this script
    # split at 0.5*RMSE_inf and reported the sub-threshold rows as "genuinely
    # usable" -- but they all sat at ~11,066 kcal/mol, which is unusable by any
    # standard. A threshold picked relative to the sentinel only tells you how
    # close to the sentinel you are. USABLE is anchored to the promotion
    # script's own MAX_ERR (100 kcal/mol), which is a real decision boundary.
    BANDS = [(0, 10, 'usable          (< 10)'),
             (10, 100, 'weak            (10-100)'),
             (100, 1000, 'very weak       (100-1k)'),
             (1000, float('inf'), 'unconstrained   (> 1k)')]
    print(f'\nNEWLY CARRYING A VALUE: {len(gained)}')
    errs = [(k, num(new[k], errcol)) for k in gained]
    errs = [(k, e) for k, e in errs if e is not None]
    for lo, hi, label in BANDS:
        sel = [k for k, e in errs if lo <= e < hi]
        print(f'    sigma kcal/mol  {label:<26}{len(sel):7d}')
        if sel and lo < 100:
            for k in sel[:5]:
                print(f'        {k}  {new[k][dgcol]:>12} +- {new[k][errcol]}')
    if errs:
        e = [x for _, x in errs]
        print(f'    sigma: min {min(e):,.1f}  median {median(e):,.1f}  '
              f'max {max(e):,.1f} kcal/mol')
        print(f'    (RMSE_inf = {RMSE_INF_KCAL:,.0f} kcal/mol = the '
              f'"unconstrained" marker)')

    # 3. agreement on shared ok rows
    both = [k for k in (o_ok & n_ok)
            if num(old[k], dgcol) is not None and num(new[k], dgcol) is not None]
    d = [abs(num(new[k], dgcol) - num(old[k], dgcol)) for k in both]
    moved = sum(1 for x in d if x > 0.01)
    print(f'\nAGREEMENT on {len(both)} rows ok in both:')
    print(f'    median |change| {median(d):.4f} kcal/mol   changed by >0.01: {moved}')

    # status census
    print('\nSTATUS CENSUS (new):')
    c = Counter(r['status'].split(':')[0] for r in new.values())
    for k, v in c.most_common():
        print(f'    {k:<46}{v:7d}')
    print()

def main():
    ap = argparse.ArgumentParser()
    MSD = (Path(__file__).resolve().parents[3] / ''
               'Biochemistry/Thermodynamics/eQuilibrator')
    E = ROOT / 'data'
    ap.add_argument('--kind', choices=('reaction', 'compound', 'both'), default='both')
    a = ap.parse_args()
    if a.kind in ('reaction', 'both'):
        report(load(MSD/'ModelSEED_Reaction_Energies.tsv', 'reaction_id'),
               load(E/'modelseed_energies.tsv', 'reaction_id'),
               'reaction_id', 'dg_prime_kcal_per_mol', 'uncertainty_kcal_per_mol',
               'REACTIONS')
    if a.kind in ('compound', 'both'):
        report(load(MSD/'ModelSEED_Compound_Energies.tsv', 'compound_id'),
               load(E/'modelseed_formation_energies.tsv', 'compound_id'),
               'compound_id', 'dgf_prime_kcal_per_mol', 'uncertainty_kcal_per_mol',
               'COMPOUNDS')

if __name__ == '__main__':
    main()
