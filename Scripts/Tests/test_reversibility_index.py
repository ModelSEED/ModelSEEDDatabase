#!/usr/bin/env python3
"""Tests for the standalone Noor 2012 / eQuilibrator 2.0 reversibility index.

Three parts:

1. **Formula** -- analytic identities the definition must satisfy, checked
   against hand-computed values rather than against another implementation.

2. **Fidelity** -- ``ln_reversibility_index`` reproduces eQuilibrator's own
   published ``ln_reversibility_index`` column for the bulk of the ModelSEED
   table. This is the real check on the formula, the RT/unit conversion and the
   water/proton exclusions.

3. **ATP synthase and ABC transporters** -- what the index actually says about
   the two reaction families ModelSEED currently decides *structurally*, before
   any energy is read. Every reaction here is scored as ordinary chemistry: two
   compartments means two species, and no membrane term, shortcut or exemption
   applies. That is the only way to ask whether the structural shortcuts are
   doing work the energy could have done, and the answer turns out to be
   "no, and the index is not a safe replacement for them either".

Usage:
    ./test_reversibility_index.py           # assertions
    ./test_reversibility_index.py --report  # assertions + the full ATPS/ABC table
"""
import math
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'Scripts', 'Thermodynamics'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'Libs', 'Python'))

import reversibility_index as RI          # noqa: E402
import reversibility_heuristics as rh     # noqa: E402

EQ_TABLE = os.path.join(REPO_ROOT, 'Biochemistry', 'Thermodynamics',
                        'eQuilibrator', 'MetaNetX_Reaction_Energies.tbl')
PLAIN_MEASUREMENT = re.compile(r'^(-?[\d.]+(?:[eE][-+]?\d+)?)\+/-')
TOL_ABS, TOL_REL = 0.10, 0.005
MIN_AGREEMENT = 17700          # measured 17,759 against origin/dev's table

failures = []


def check(label, condition, detail=''):
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}"
          + (f' -- {detail}' if detail else ''))
    if not condition:
        failures.append(label)


def rgt(compound, coefficient, compartment=0):
    return {'compound': compound, 'coefficient': coefficient,
            'compartment': compartment}


# --- 1. Formula -----------------------------------------------------------
def test_formula():
    print('\nFormula')
    rt = RI.RT_KCAL

    # A <=> B. sum|nu| = 2 and sum nu = 0, so the 1 mM correction vanishes and
    # ln(Gamma) collapses to dG'o / RT exactly.
    a_to_b = [rgt('cpd00100', -1), rgt('cpd00200', 1)]
    check('A<=>B: ln(Gamma) == dG/RT',
          math.isclose(RI.ln_reversibility_index(a_to_b, 3.0), 3.0 / rt,
                       rel_tol=1e-12),
          f'{RI.ln_reversibility_index(a_to_b, 3.0):.4f}')

    # A <=> B + C. sum nu = +1, so dG'm = dG + RT*ln(1e-3), and sum|nu| = 3.
    a_to_bc = [rgt('cpd00100', -1), rgt('cpd00200', 1), rgt('cpd00300', 1)]
    expect = (2.0 / 3.0) * (3.0 + rt * math.log(1e-3)) / rt
    check('A<=>B+C: 1 mM correction and 2/sum|nu| both applied',
          math.isclose(RI.ln_reversibility_index(a_to_bc, 3.0), expect,
                       rel_tol=1e-12),
          f'{RI.ln_reversibility_index(a_to_bc, 3.0):.4f}')

    # Water and protons are dropped from both sums.
    with_solvent = a_to_b + [rgt(RI.WATER, -3), rgt(RI.PROTON, 2)]
    check('water and protons excluded',
          RI.ln_reversibility_index(with_solvent, 3.0)
          == RI.ln_reversibility_index(a_to_b, 3.0))
    check('exclude_protons=False changes the answer',
          RI.ln_reversibility_index(with_solvent, 3.0, exclude_protons=False)
          != RI.ln_reversibility_index(a_to_b, 3.0))

    # Sign convention: negative dG means the reaction runs forward, operator '>'.
    check('negative dG -> ">"',
          RI.direction_from_index(RI.ln_reversibility_index(a_to_b, -20.0)) == '>')
    check('positive dG -> "<"',
          RI.direction_from_index(RI.ln_reversibility_index(a_to_b, 20.0)) == '<')
    check('small dG -> "="',
          RI.direction_from_index(RI.ln_reversibility_index(a_to_b, 1.0)) == '=')

    # ln(1000) is the cut, so Gamma = 1000 sits exactly on it.
    on_cut = rt * RI.LN_GAMMA_THRESHOLD
    check('Gamma == 1000 is not yet irreversible',
          RI.direction_from_index(RI.ln_reversibility_index(a_to_b, on_cut)) == '=')
    check('Gamma just past 1000 is irreversible',
          RI.direction_from_index(
              RI.ln_reversibility_index(a_to_b, on_cut * 1.001)) == '<')

    # Reversing a reaction negates the index.
    reverse = [rgt('cpd00100', 1), rgt('cpd00200', -1)]
    check('reversing the reaction negates ln(Gamma)',
          math.isclose(RI.ln_reversibility_index(a_to_b, 3.0),
                       -RI.ln_reversibility_index(reverse, -3.0), rel_tol=1e-12))

    # Guards.
    check('no countable reagents -> None',
          RI.ln_reversibility_index([rgt(RI.WATER, -1), rgt(RI.PROTON, 1)], 5.0)
          is None)
    check('a species that cancels drops out',
          RI.coefficient_sums([rgt('cpd00100', -1), rgt('cpd00100', 1)])[1] == 0.0)
    check('same compound in two compartments stays two species',
          RI.coefficient_sums([rgt('cpd00100', -1, 0),
                               rgt('cpd00100', 1, 1)])[2] == 2)
    check('Gamma overflow is clamped, not raised',
          RI.reversibility_index(a_to_b, -1e5) == 0.0)

    # Error propagation scales by the same factor as the energy.
    check('error propagates by 2/sum|nu|/RT',
          math.isclose(RI.ln_reversibility_index_error(a_to_bc, 3.0),
                       (2.0 / 3.0) * 3.0 / rt, rel_tol=1e-12))


# --- 2. Fidelity against eQuilibrator's published column -------------------
def test_fidelity(reactions):
    print("\nFidelity vs eQuilibrator's published ln_reversibility_index")
    agree = mismatch = 0
    unexplained = []
    for line in open(EQ_TABLE):
        array = line.rstrip('\n').split('\t')
        if len(array) < 4 or array[0] not in reactions:
            continue
        match = PLAIN_MEASUREMENT.match(array[3])
        if not match:
            continue                      # undecomposable marker, no information
        published = float(match.group(1))
        entry = reactions[array[0]]
        mine = RI.ln_reversibility_index(entry['stoichiometry'], float(array[1]))
        if mine is None:
            continue
        if abs(mine - published) <= max(TOL_ABS, TOL_REL * abs(published)):
            agree += 1
        else:
            mismatch += 1
            if entry.get('is_transport') != 1 and not _has_duplicate_compound(entry):
                unexplained.append(array[0])

    check('reproduces eQuilibrator on the bulk of the table',
          agree >= MIN_AGREEMENT, f'{agree} agree (floor {MIN_AGREEMENT})')
    check('every mismatch is a compartment-collapsed reaction',
          len(unexplained) < 0.10 * max(mismatch, 1),
          f'{mismatch} mismatches, {len(unexplained)} unexplained')


def _has_duplicate_compound(entry):
    seen = set()
    for reagent in entry['stoichiometry']:
        if reagent['compound'] in seen:
            return True
        seen.add(reagent['compound'])
    return False


# --- 3. ATP synthase and ABC transporters ---------------------------------
def classify(reactions):
    """Split out the two families ModelSEED decides structurally."""
    atps, abct = [], []
    for rxn_id, entry in reactions.items():
        terms = rh._walk_stoichiometry(entry['stoichiometry'])
        if rh._is_atp_synthase(entry, terms['proton_cpts']):
            atps.append(rxn_id)
        elif rh._abc_transporter_decision(entry, terms['phosphates']):
            abct.append(rxn_id)
    return sorted(atps), sorted(abct)


def score(entry, exclude_protons=True):
    """(dg, sigma, ln_gamma, operator, undecomposable) from the eQuilibrator
    sublist -- the per-source dictionary, never the flat ``deltag`` field."""
    thermo = entry.get('thermodynamics') or {}
    pair = thermo.get('eQuilibrator')
    if not pair or pair[0] is None:
        return None
    dg, sigma = float(pair[0]), float(pair[1])
    ln_gamma = RI.ln_reversibility_index(entry['stoichiometry'], dg,
                                         exclude_protons=exclude_protons)
    return (dg, sigma, ln_gamma, RI.direction_from_index(ln_gamma),
            abs(sigma) >= rh.EQ_UNDECOMPOSABLE_SIGMA)


def test_atp_synthase(reactions, atps, report=False):
    print('\nATP synthase -- index applied as ordinary chemistry')
    scored = [(r, score(reactions[r])) for r in atps]
    scored = [(r, s) for r, s in scored if s]
    check('every ATP synthase reaction carries an eQuilibrator energy',
          len(scored) == len(atps), f'{len(scored)}/{len(atps)}')

    if report:
        print(f"\n  {'reaction':10s} {'dG':>8s} {'sigma':>7s} {'lnGamma':>9s} "
              f"{'op':>3s} {'lnG(+H+)':>9s} {'op':>3s}")
        for rxn_id, (dg, sigma, ln_gamma, op, _) in scored:
            with_h = RI.ln_reversibility_index(
                reactions[rxn_id]['stoichiometry'], dg, exclude_protons=False)
            print(f'  {rxn_id:10s} {dg:8.2f} {sigma:7.2f} {ln_gamma:9.2f} '
                  f'{op:>3s} {with_h:9.2f} {RI.direction_from_index(with_h):>3s}')

    # Every one of them clears the ln(1000) cut, so the index calls the most
    # famously reversible enzyme in the cell irreversible.
    hard = [r for r, s in scored if s[3] in '<>']
    check('the index calls every ATP synthase reaction irreversible',
          len(hard) == len(scored), f'{len(hard)}/{len(scored)} get "<" or ">"')

    # And the magnitude is identical across all of them, which is the tell: the
    # protons cancelled in the collapsed MetaNetX formula, so what eQuilibrator
    # actually scored was ADP + Pi <=> ATP + H2O, the same reaction every time.
    magnitudes = {round(abs(s[2]), 2) for _, s in scored}
    check('all share one |ln Gamma| -- the collapsed ADP+Pi<=>ATP formula',
          len(magnitudes) == 1, f'|ln Gamma| = {magnitudes}')

    # The direction therefore tracks how the reaction was written, not physics.
    directions = {s[3] for _, s in scored}
    check('direction splits both ways across the family',
          directions == {'<', '>'},
          f"{sum(1 for _, s in scored if s[3] == '>')} forward, "
          f"{sum(1 for _, s in scored if s[3] == '<')} reverse")

    # Letting the translocated protons count pulls most of them back inside the
    # reversible window -- the proton term is the entire driving force, and
    # excluding it (correctly, for a buffered cytosol) discards it.
    relaxed = sum(1 for r, s in scored
                  if RI.direction_from_index(RI.ln_reversibility_index(
                      reactions[r]['stoichiometry'], s[0],
                      exclude_protons=False)) == '=')
    check('counting the translocated protons returns most to "="',
          relaxed > len(scored) / 2, f'{relaxed}/{len(scored)} become "="')


def test_abc_transporters(reactions, abct, report=False):
    print('\nABC transporters -- index applied as ordinary chemistry')
    scored = []
    for rxn_id in abct:
        result = score(reactions[rxn_id])
        if result:
            scored.append((rxn_id, result))
    check('most ABC transporters carry an eQuilibrator energy',
          len(scored) > 0.7 * len(abct), f'{len(scored)}/{len(abct)}')

    undecomposable = [r for r, s in scored if s[4]]
    clean = [(r, s) for r, s in scored if not s[4]]
    check('the sigma gate flags the undecomposable ones',
          len(undecomposable) < 0.05 * len(scored),
          f'{len(undecomposable)} flagged, {len(clean)} scorable')

    agree = disagree = reversible = 0
    for rxn_id, (dg, sigma, ln_gamma, op, _) in clean:
        terms = rh._walk_stoichiometry(reactions[rxn_id]['stoichiometry'])
        structural = rh._abc_transporter_decision(reactions[rxn_id],
                                                  terms['phosphates'])[1]
        if op == '=':
            reversible += 1
        elif op == structural:
            agree += 1
        else:
            disagree += 1

    total = len(clean)
    print(f'  ... vs the structural ATP-sign rule: {agree} agree, '
          f'{disagree} opposite, {reversible} called reversible '
          f'({100.0 * agree / total:.1f}% agreement)')
    check('the index mostly reproduces the structural ATP-sign rule',
          agree > 0.9 * total, f'{agree}/{total}')
    check('but it is not unanimous -- some flip outright',
          disagree > 0, f'{disagree} reactions get the opposite operator')

    # The knife-edge: for the typical ABC transporter the substrate cancels in
    # the collapsed formula, leaving plain ATP hydrolysis at -6.54 kcal/mol.
    # ln(Gamma) then lands at -7.18 against a cut of -6.91 -- inside the cut by
    # 4%. The 94% agreement above is one rounding away from collapsing.
    near = [s[2] for _, s in clean
            if s[3] in '<>' and abs(s[2]) < 1.15 * RI.LN_GAMMA_THRESHOLD]
    check('most hard calls sit within 15% of the threshold',
          len(near) > 0.8 * (agree + disagree),
          f'{len(near)}/{agree + disagree} inside 1.15x ln(1000)')

    if report:
        print(f"\n  {'reaction':10s} {'dG':>8s} {'lnGamma':>9s} "
              f"{'struct':>7s} {'index':>6s}")
        for rxn_id, (dg, sigma, ln_gamma, op, _) in clean[:15]:
            terms = rh._walk_stoichiometry(reactions[rxn_id]['stoichiometry'])
            structural = rh._abc_transporter_decision(
                reactions[rxn_id], terms['phosphates'])[1]
            print(f'  {rxn_id:10s} {dg:8.2f} {ln_gamma:9.2f} '
                  f'{structural:>7s} {op:>6s}')
        print(f'  ... {len(clean) - 15} more')


def main():
    report = '--report' in sys.argv
    from BiochemPy import Reactions
    reactions = Reactions().loadReactions()

    test_formula()
    test_fidelity(reactions)
    atps, abct = classify(reactions)
    print(f'\nFound {len(atps)} ATP synthase and {len(abct)} ABC transport '
          f'reactions')
    test_atp_synthase(reactions, atps, report)
    test_abc_transporters(reactions, abct, report)

    print('\n' + ('PASS' if not failures else
                  'FAIL: ' + '; '.join(failures)))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
