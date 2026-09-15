#!/usr/bin/env python
"""Set the canonical `reversibility` field from the GRADED RECOMMENDATION.

Written 2026-09-15, replacing the fixed source precedence in
Estimate_Reaction_Reversibility.py for the canonical field.

WHY. The old rule read eQuilibrator's energy first, whatever the grading had
decided. On the 3,386 reactions where eQuilibrator disclaims (sigma >= 2500,
its undecomposable marker) that meant the canonical field was computed from a
source the grading had already vetoed -- so a reaction could be
`silver, dGPredictor, forward` in thermo-evidence and `?` in reversibility at
the same time. 1,230 of those carry a committed direction from their graded
source that the canonical field contradicted.

THE RULE. One mechanism now produces both:

    graded          -> the operator of the reaction's best_source
    best_source is  -> computed from the experimental dG and its sd, via the
      openTECR         same per-source operator cascade (it has no entry in
                       `thermodynamics`, carrying a measurement not a
                       prediction)
    ungraded        -> '?'   no source has a usable energy, so there is
                       nothing to base a direction on

This deliberately supersedes the 2020 canonical semantics. The effect on
downstream models is to be quantified separately.
"""
import sys, os, csv, json
sys.path.append('../../Libs/Python/')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from BiochemPy import Reactions
import _thermo_helpers as th

GRADES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', '..', 'Biochemistry', 'Thermodynamics',
                      'SourceGrading', 'results', 'thermo_grades')
KEY = {'EQ': 'eQuilibrator', 'DGP': 'dGPredictor',
       'GC': 'Group contribution', 'openTECR': 'openTECR'}


def main():
    dry = '--dry-run' in sys.argv
    with open(os.path.join(GRADES, 'reaction_grades.tsv')) as fh:
        rg = {r['rxn']: r for r in csv.DictReader(fh, delimiter='\t')}
    meas = {}
    with open(os.path.join(GRADES, 'source_grades.tsv')) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['source'] == 'openTECR' and r['dg']:
                meas[r['rxn']] = (float(r['dg']),
                                  float(r['sigma']) if r['sigma'] else 0.0)

    helper = Reactions()
    rxns = helper.loadReactions()
    from collections import Counter
    n = Counter()
    for rid, entry in rxns.items():
        g = rg.get(rid, {})
        grade = g.get('best_grade')
        if not grade:
            op = '?'; n['ungraded -> ?'] += 1
        else:
            label = KEY.get(g.get('best_source', ''), '')
            if label == 'openTECR':
                dg, sd = meas.get(rid, (None, None))
                op = (th._per_source_operator(entry, dg, sd, 'openTECR')
                      if dg is not None else '?')
                n['openTECR (from measurement)'] += 1
            else:
                trip = (entry.get('thermodynamics') or {}).get(label)
                op = trip[2] if isinstance(trip, list) and len(trip) > 2 else '?'
                n[f'{label}'] += 1
        if entry.get('reversibility') != op:
            n['CHANGED'] += 1
        entry['reversibility'] = op

    for k, v in sorted(n.items(), key=lambda t: -t[1]):
        print(f'  {k:<34} {v:>7,}')
    print(f'  {"final distribution":<34} '
          f'{dict(Counter(e.get("reversibility") for e in rxns.values()))}')
    if dry:
        print('  DRY RUN -- nothing written.')
        return
    helper.saveReactions(rxns)


if __name__ == '__main__':
    main()
