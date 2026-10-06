#!/usr/bin/env python3
"""Landmark regression test for the Group-Contribution Convention A pipeline.

Guards against three classes of silent regression:

1. **Anchor drift** — the 9 Convention A load-bearing constants (H+, H2O,
   NH4+, CO2, HCO3, H2O2, H2, O2, H2S) are checked exactly. If any anchor
   changes value, someone edited the anchor table or the injection step
   silently stopped running.
2. **Landmark drift** — a dozen central-metabolism compounds (ATP, NAD,
   glucose, ...) have hardcoded expected ΔG values. These must match to
   ±0.05 kcal/mol. Bigger drift means the pipeline itself changed
   (different resolver, mean vs lowest, new MolAnalysis inputs, etc.).
3. **Convention flip** — water at -37.54 kcal/mol would indicate the DB
   accidentally went back to Convention B (Alberty) instead of Convention A
   (Chris). Explicit assertion catches this.

Also spot-checks two landmark reactions (rxn00001, rxn00566) to catch
regressions where compound values are right but the aggregation is broken.

Exit code 0 on all-pass, 1 on any failure. Suitable for CI.
"""
import json
import os
import sys
from glob import glob

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
BIOCHEM_DIR = os.path.join(THIS_DIR, '..', '..', 'Biochemistry')

LABEL = 'Group contribution'
ANCHOR_TOL = 0.01     # anchors are hardcoded, must match exactly (float rounding)
LANDMARK_TOL = 0.05   # landmarks are pipeline output, tiny rounding slack
REACTION_TOL = 0.5    # reactions accumulate rounding across many compounds
SENTINEL = 10000000.0

# The 9 Convention A anchors — cf.
# Scripts/Thermodynamics/Update_Compound_GroupContribution_Energies.py
# CONVENTION_A_ANCHORS. Values must exactly match what's stamped.
ANCHORS = [
    # (cpd_id, name,   expected_dg, expected_dge)
    ('cpd00067', 'H+',    -9.5,      0.0),
    ('cpd00001', 'H2O',  -56.687,    0.0),
    ('cpd00013', 'NH4+', -18.97,     0.0),
    ('cpd00011', 'CO2',  -92.26,     0.0),
    ('cpd00242', 'HCO3', -140.26,    0.0),
    ('cpd00025', 'H2O2', -32.05,     0.0),
    ('cpd11640', 'H2',     4.2065,   0.0),
    ('cpd00007', 'O2',     3.9197,   0.0),
    ('cpd00239', 'H2S',   -6.66,     0.0),
]

# Central-metabolism landmarks. Values captured 2026-08-07 after the
# Phase-3 rebuild (mean-of-aliases resolver over regenerated MolAnalysis
# tables with the 30s cycle-timeout + sentinel-on-timeout patch).
# Any change here means the pipeline itself moved.
LANDMARKS = [
    # (cpd_id, name,       expected_dg,  expected_dge)
    ('cpd00002', 'ATP',        -673.85,   6.09),
    ('cpd00003', 'NAD',        -529.59,   8.71),
    ('cpd00004', 'NADH',       -524.32,   8.54),
    ('cpd00005', 'NADPH',      -736.82,   8.52),
    ('cpd00006', 'NADP',       -742.09,   8.69),
    ('cpd00008', 'ADP',        -465.85,   6.07),
    ('cpd00009', 'Pi',         -261.97,   1.0),
    ('cpd00010', 'CoA',        -751.99,   6.9),
    ('cpd00012', 'PPi',        -480.93,   1.0),
    ('cpd00020', 'Pyruvate',   -112.69,   0.6),
    ('cpd00023', 'L-Glutamate',-164.13,   0.75),
    ('cpd00027', 'D-Glucose',  -218.28,   3.09),
]

# Reaction spot-checks — GC value should aggregate consistently from
# the per-compound values above.
REACTION_LANDMARKS = [
    # (rxn_id, name,           expected_dg, expected_dge)
    ('rxn00001', 'PPi hydrolysis',         4.18,   2.24),
    ('rxn00566', 'Cysteine desulfhydration', -11.68, 5.76),
]

# Convention A sentinel value that must NOT appear for water. If it does,
# the DB slipped back to Alberty/Convention B and every ΔfG is off by
# n_H × 9.539 kcal/mol for a given compound.
CONVENTION_B_WATER = -37.6   # rough — actual Alberty value; +/- 1 kcal/mol tolerance


def load_compounds():
    out = {}
    for path in sorted(glob(os.path.join(BIOCHEM_DIR, 'compound_*.json'))):
        with open(path) as fh:
            for c in json.load(fh):
                out[c['id']] = c
    return out


def load_reactions():
    out = {}
    for path in sorted(glob(os.path.join(BIOCHEM_DIR, 'reaction_*.json'))):
        with open(path) as fh:
            for r in json.load(fh):
                out[r['id']] = r
    return out


def get_gc(entry):
    """Return (dg, dge) or None if no numeric GC entry."""
    thermo = entry.get('thermodynamics')
    if not isinstance(thermo, dict):
        return None
    gc = thermo.get(LABEL)
    if not isinstance(gc, list) or len(gc) < 2:
        return None
    try:
        dg = float(gc[0])
        dge = float(gc[1])
    except (TypeError, ValueError):
        return None
    return (dg, dge)


def check_anchors(compounds):
    failures = []
    for cpd_id, name, exp_dg, exp_dge in ANCHORS:
        entry = compounds.get(cpd_id)
        if entry is None:
            failures.append(f'  {cpd_id} ({name}): MISSING from compound_*.json')
            continue
        gc = get_gc(entry)
        if gc is None:
            failures.append(f'  {cpd_id} ({name}): no numeric GC entry (should be {exp_dg})')
            continue
        dg, dge = gc
        if abs(dg - exp_dg) > ANCHOR_TOL or abs(dge - exp_dge) > ANCHOR_TOL:
            failures.append(
                f'  {cpd_id} ({name}): expected [{exp_dg}, {exp_dge}], '
                f'got [{dg}, {dge}]')
    return failures


def check_landmarks(compounds):
    failures = []
    for cpd_id, name, exp_dg, exp_dge in LANDMARKS:
        entry = compounds.get(cpd_id)
        if entry is None:
            failures.append(f'  {cpd_id} ({name}): MISSING')
            continue
        gc = get_gc(entry)
        if gc is None:
            failures.append(f'  {cpd_id} ({name}): no numeric GC (expected {exp_dg})')
            continue
        dg, dge = gc
        if abs(dg - exp_dg) > LANDMARK_TOL:
            failures.append(
                f'  {cpd_id} ({name}): dg drifted — expected {exp_dg}, got {dg} '
                f'(|diff|={abs(dg-exp_dg):.3f} > {LANDMARK_TOL})')
        if abs(dge - exp_dge) > LANDMARK_TOL:
            failures.append(
                f'  {cpd_id} ({name}): dge drifted — expected {exp_dge}, got {dge}')
    return failures


def check_reactions(reactions):
    failures = []
    for rxn_id, name, exp_dg, exp_dge in REACTION_LANDMARKS:
        entry = reactions.get(rxn_id)
        if entry is None:
            failures.append(f'  {rxn_id} ({name}): MISSING')
            continue
        gc = get_gc(entry)
        if gc is None:
            failures.append(f'  {rxn_id} ({name}): no numeric GC (expected {exp_dg})')
            continue
        dg, dge = gc
        if abs(dg - exp_dg) > REACTION_TOL:
            failures.append(
                f'  {rxn_id} ({name}): dg drifted — expected {exp_dg}, got {dg}')
    return failures


def check_convention_consistency(compounds):
    """Water is the canary. Convention A says -56.687; Convention B (Alberty)
    would say ~-37.6. If water lands closer to -37.6 than -56.687, the DB
    silently flipped conventions and every H-containing compound is off."""
    failures = []
    water = compounds.get('cpd00001')
    gc = get_gc(water) if water else None
    if gc is None:
        failures.append('  cpd00001 (H2O): no GC entry — cannot verify convention')
        return failures
    dg = gc[0]
    dist_A = abs(dg - (-56.687))
    dist_B = abs(dg - CONVENTION_B_WATER)
    if dist_B < dist_A:
        failures.append(
            f'  cpd00001 (H2O): dg={dg} — CLOSER TO CONVENTION B (Alberty, {CONVENTION_B_WATER}) '
            f'than to Convention A (-56.687). Convention has silently flipped; '
            f'every H-containing compound is now off by n_H × 9.539 kcal/mol.')
    return failures


def main():
    print('Loading compounds and reactions...')
    compounds = load_compounds()
    reactions = load_reactions()
    print(f'  Loaded {len(compounds):,} compounds, {len(reactions):,} reactions\n')

    all_failures = []

    print('── Convention A anchors (must match exactly to ±0.01 kcal/mol) ──')
    f = check_anchors(compounds)
    if f:
        print('  FAIL')
        for line in f:
            print(line)
        all_failures.extend(f)
    else:
        print(f'  PASS — all {len(ANCHORS)} anchors match')

    print('\n── Central-metabolism landmarks (must match to ±0.05 kcal/mol) ──')
    f = check_landmarks(compounds)
    if f:
        print('  FAIL')
        for line in f:
            print(line)
        all_failures.extend(f)
    else:
        print(f'  PASS — all {len(LANDMARKS)} landmarks match')

    print('\n── Reaction spot-checks (must match to ±0.5 kcal/mol) ──')
    f = check_reactions(reactions)
    if f:
        print('  FAIL')
        for line in f:
            print(line)
        all_failures.extend(f)
    else:
        print(f'  PASS — all {len(REACTION_LANDMARKS)} reaction landmarks match')

    print('\n── Convention A vs B consistency (water canary) ──')
    f = check_convention_consistency(compounds)
    if f:
        print('  FAIL')
        for line in f:
            print(line)
        all_failures.extend(f)
    else:
        print('  PASS — water at Convention A (H included)')

    print()
    if all_failures:
        print(f'RESULT: FAIL — {len(all_failures)} issue(s)')
        sys.exit(1)
    print('RESULT: PASS')
    sys.exit(0)


if __name__ == '__main__':
    main()
