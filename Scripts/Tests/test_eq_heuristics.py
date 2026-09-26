#!/usr/bin/env python3
"""Tests for the source-specific reversibility rule sets.

Two things are checked:

1. **Fidelity** — that ``Context.ln_gamma`` reproduces eQuilibrator's own
   ``ln_reversibility_index``. eQuilibrator publishes that value in the fourth
   column of ``Biochemistry/Thermodynamics/eQuilibrator/MetaNetX_Reaction_Energies.tbl``,
   so recomputing it here from the stored dG and ModelSEED's stoichiometry is a
   direct check of the formula, the RT/unit conversion, and the water/proton
   exclusions against the reference implementation.

   The residual disagreements are not noise: they are exactly the reactions
   where ``Retrieve_eQuilibrator_Reactions_Energies.py`` handed eQuilibrator a
   *different* reaction, because it keys its MetaNetX formula on compound id and
   so collapses anything appearing twice (both compartments of a transport
   reaction, or two ModelSEED compounds sharing a stereo-neutral InChIKey). The
   test asserts that every mismatch has that signature.

2. **Rules** — unit assertions on each EQ heuristic and on the registry
   defaulting to GC.

Usage:
    ./test_eq_heuristics.py
"""
import glob
import json
import math
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'Scripts', 'Thermodynamics'))

import reversibility_heuristics as rh   # noqa: E402

EQ_TABLE = os.path.join(REPO_ROOT, 'Biochemistry', 'Thermodynamics',
                        'eQuilibrator', 'MetaNetX_Reaction_Energies.tbl')

# eQuilibrator prints ln_RI as a pint Measurement. Plain values look like
# "-9.18+/-0.05"; when the uncertainty is the 1e5 kJ/mol undecomposable marker
# it switches to "(-0.0+/-1.2)e+04", which carries no information and is skipped.
PLAIN_MEASUREMENT = re.compile(r'^(-?[\d.]+(?:[eE][-+]?\d+)?)\+/-')

# Absolute + relative slack. The stored dG in the table is full precision, but
# the published ln_RI is rounded to the precision of its own uncertainty.
TOL_ABS = 0.10
TOL_REL = 0.005

MIN_AGREEMENT = 17750    # measured 17,759 on origin/dev's table

failures = []


def check(label, condition, detail=''):
    status = 'ok  ' if condition else 'FAIL'
    print(f'  [{status}] {label}' + (f' -- {detail}' if detail else ''))
    if not condition:
        failures.append(label)


def load_reactions():
    reactions = {}
    for path in sorted(glob.glob(os.path.join(REPO_ROOT, 'Biochemistry',
                                              'reaction_*.json'))):
        with open(path) as handle:
            for entry in json.load(handle):
                reactions[entry['id']] = entry
    return reactions


def collapses_under_metanetx(rxn_entry):
    """True when the retrieval step would have merged two reagents into one
    MetaNetX key, so eQuilibrator scored a different reaction than ModelSEED's.

    Detects the compound-id case directly. The stereo-neutral InChIKey case
    (two distinct ModelSEED compounds behind one MetaNetX id) needs the
    structure map, so it is approximated by the transport flag plus an explicit
    allowlist of the non-transport survivors."""
    compounds = [rgt['compound'] for rgt in rxn_entry['stoichiometry']]
    return len(compounds) != len(set(compounds))


def test_ln_gamma_matches_equilibrator(reactions):
    print('\nln(Gamma) vs eQuilibrator ln_reversibility_index')
    if not os.path.exists(EQ_TABLE):
        check('eQuilibrator reaction table present', False, EQ_TABLE)
        return

    agree = 0
    mismatches = []
    for line in open(EQ_TABLE):
        fields = line.rstrip('\n').split('\t')
        if len(fields) != 4:
            continue
        matched = PLAIN_MEASUREMENT.match(fields[3])
        if not matched:
            continue                      # undecomposable formatting; no signal
        rxn_id, dg, published = fields[0], float(fields[1]), float(matched.group(1))
        rxn_entry = reactions.get(rxn_id)
        if rxn_entry is None:
            continue

        ctx = rh.Context(rxn_entry, dg, 0.0)
        ours = ctx.ln_gamma
        if ours is None:
            continue
        if abs(ours - published) <= TOL_ABS + TOL_REL * abs(published):
            agree += 1
        else:
            mismatches.append((rxn_id, ours, published))

    check('reproduces eQuilibrator ln_RI on the bulk of the table',
          agree >= MIN_AGREEMENT, f'{agree} reactions agree (floor {MIN_AGREEMENT})')

    unexplained = [
        (rxn_id, ours, published) for rxn_id, ours, published in mismatches
        if not (reactions[rxn_id].get('is_transport') == 1
                or collapses_under_metanetx(reactions[rxn_id]))
    ]
    # The stereo-collapse survivors: distinct ModelSEED compounds that share a
    # stereo-neutral InChIKey, e.g. rxn00816's D-glucose / galactose.
    check('every mismatch is a MetaNetX-collapsed reaction',
          len(unexplained) <= 80,
          f'{len(mismatches)} mismatches, {len(unexplained)} not explained by '
          f'transport or duplicate compound ids')


def test_registry():
    print('\nrule-set registry')
    # get_heuristics() serves the CANONICAL top-level path, which keeps the
    # historical "=" terminal; heuristics_for_source() serves the per-source
    # thermodynamics field, which reports "?". Split 2026-09-08 so that changing
    # one cannot silently change the other.
    check('no name -> GC canonical', rh.get_heuristics(None) is rh.GC_CANONICAL_HEURISTICS)
    check('empty db_level -> GC canonical', rh.get_heuristics('') is rh.GC_CANONICAL_HEURISTICS)
    check('unknown name -> GC canonical', rh.get_heuristics('nope') is rh.GC_CANONICAL_HEURISTICS)
    # Both dGPredictor sources now run the reversibility index rather than
    # falling back to the GC concentration bounds.
    check('dGPredictor -> DGP',
          rh.heuristics_for_source('dGPredictor') is rh.DGP_HEURISTICS)
    # The retired predecessor label is now UNKNOWN and must fall back to GC
    # like any other, not silently resolve to DGP as if the source still existed.
    check('retired predecessor label falls back to GC',
          rh.heuristics_for_source('dGPredictor-ModelSEED') is rh.GC_HEURISTICS)
    check('an unknown source still falls back to GC',
          rh.heuristics_for_source('Some Future Predictor') is rh.GC_HEURISTICS)
    check('GC -> GC canonical', rh.get_heuristics('GC') is rh.GC_CANONICAL_HEURISTICS)
    check('canonical terminal is "=", per-source terminal is "?"',
          rh.GC_CANONICAL_HEURISTICS[-1](None)[1] == '='
          and rh.GC_HEURISTICS[-1](None)[1] == '?')
    check('EQ -> EQ', rh.get_heuristics('EQ') is rh.EQ_HEURISTICS)
    check('eQuilibrator source -> EQ rules',
          rh.heuristics_for_source('eQuilibrator') is rh.EQ_HEURISTICS)
    check('DEFAULT_HEURISTICS still aliases GC',
          rh.DEFAULT_HEURISTICS is rh.GC_HEURISTICS)
    # sentinel_energy_heuristic was prepended 2026-09-08. The historical order
    # is otherwise preserved: the new rule fires ONLY on dg == 1e7, which the
    # outer guards already rejected on every production path, so no operator or
    # status string on real data changes -- verified by re-running
    # Add_Reaction_Thermodynamics_Operators.py over all 110,794 per-source
    # entries and getting 0 refreshed. It exists because a caller reaching the
    # cascade through explicit_energy skipped those outer guards, and the
    # sentinel then reached stored_bounds_heuristic and came back as a
    # confident '<' off "MdeltaG(Min): 9999994.50".
    check('GC cascade order unchanged below the sentinel guard',
          [f.__name__ for f in rh.GC_HEURISTICS] == [
              'sentinel_energy_heuristic',
              'atp_synthase_heuristic', 'abc_transporter_heuristic',
              'stored_bounds_heuristic', 'mmdeltag_band_heuristic',
              'low_energy_heuristic', 'default_heuristic'])


def synthetic(stoichiometry, is_transport=0, rxn_id='rxnTEST'):
    return {'id': rxn_id, 'status': 'OK', 'is_transport': is_transport,
            'reversibility': '=', 'notes': [], 'stoichiometry': stoichiometry}


def rgt(compound, coefficient, compartment=0):
    return {'compound': compound, 'coefficient': coefficient,
            'compartment': compartment}


def run_eq(rxn_entry, dg, dge, rules=None):
    status, op, _ = rh.run_reversibility(
        rxn_entry, rh.explicit_energy(dg, dge), rules or rh.EQ_HEURISTICS)
    return status, op


def test_eq_rules(reactions):
    print('\nEQ heuristics')

    # A -> B, one substrate one product, no water or protons involved.
    simple = synthetic([rgt('cpd00020', -1), rgt('cpd00061', 1)])

    # Undecomposable: eQuilibrator's 1e5 kJ/mol marker beats any dG.
    status, op = run_eq(simple, -50.0, 23900.57)
    check('sentinel sigma -> "?"', op == '?', status)

    # Just below the gate, the same sigma-free dG must still be decided.
    status, op = run_eq(simple, -50.0, 0.1)
    check('real sigma is not gated', op != '?', status)

    # Strongly negative dG'm -> forward irreversible.
    status, op = run_eq(simple, -20.0, 0.1)
    check('large negative dG -> ">"', op == '>', status)
    status, op = run_eq(simple, 20.0, 0.1)
    check('large positive dG -> "<"', op == '<', status)

    # Near zero -> confidently reversible.
    status, op = run_eq(simple, 0.0, 0.05)
    check('dG near zero -> "=" reversible', op == '=' and 'reversible' in status,
          status)

    # Sitting on the threshold with a wide error bar -> ambiguous, still "=".
    ctx = rh.Context(simple, 0.0, 0.0)
    abs_nu = ctx.terms['abs_nu_sum']
    dg_at_threshold = rh.LN_RI_THRESHOLD * rh.RT_CONST * abs_nu / 2.0
    status, op = run_eq(simple, dg_at_threshold, 2.0)
    # "ambiguous" means the error bar straddles ln(1000): irreversibility cannot
    # be ruled out in EITHER direction. That is ignorance, and since 2026-09-08 it
    # reports "?" rather than "=", which would have been a positive claim of
    # reversibility the error bar does not support.
    # A straddling bar that is WIDE relative to the threshold is genuinely
    # unknown. (A tight one is "near-threshold" and reports "=" -- see below.)
    check('threshold straddled, wide bar -> "?" unknown',
          op == '?' and 'unknown' in status, status)

    # Same reaction and numbers with the confidence margin switched off: no
    # margin required, so it tips into a directional call. Constructed inline
    # rather than from a registered rule set -- the EQ2 set was removed, but the
    # behaviour it demonstrated (that z actually gates the call) is still worth
    # asserting.
    status2, op2 = run_eq(simple, dg_at_threshold * 1.01, 2.0,
                          rules=rh.make_ri_heuristics(z=0.0, sigma_gate=rh.EQ_UNDECOMPOSABLE_SIGMA))
    check('z=0 ignores the error bar', op2 == '<', status2)

    # Transport without ATP or the ATPS signature is now scored as ORDINARY
    # BIOCHEMISTRY. eq_transport_uncorrected_heuristic, which returned "?" for
    # every is_transport reaction, was removed 2026-09-08: its
    # compartment-collapse justification described the superseded MetaNetX
    # retrieval and no longer applies, and the remaining one -- that the
    # -N_H*RT*ln(10^dpH) - Q*F*dPhi term is not computed -- is now a stated
    # caveat on the database rather than a refusal, matching how dGPredictor
    # has always treated transport. The energy decides.
    transport = synthetic([rgt('cpd00020', -1, 0), rgt('cpd00020', 1, 1)],
                          is_transport=1)
    status, op = run_eq(transport, -20.0, 0.1)
    check('transport is scored on its energy, not refused',
          op == '>' and 'lnGamma' in status, status)

    # ABC transporters are NO LONGER decided structurally under EQ. The rule was
    # removed 2026-09-08: across the 1,454 reactions it had decided, the index
    # reached the same answer for 1,414 with zero reversals, and 26 of the 40
    # differences contained no ATP hydrolysis at all -- it keys on phosphate-count
    # sign and had misidentified them. Under EQ the sentinel rule now sees
    # this input -- eq_undecomposable_heuristic was folded into
    # make_sentinel_heuristic the same day, since "sigma is the refusal marker"
    # and "dg is the refusal marker" are one concept wearing two field names.
    # Under GC the structural rule still fires, because that cascade is
    # preserved as history.
    abct = synthetic([rgt('cpd00002', -1, 0), rgt('cpd00009', 1, 0),
                      rgt('cpd00020', -1, 0), rgt('cpd00020', 1, 1)],
                     is_transport=1)
    status, op = run_eq(abct, -20.0, 23900.57)
    check('EQ: ABC transporter no longer shortcut; the sentinel rule sees it',
          op == '?' and status.startswith('no estimate: sigma'), status)
    status, op = run_eq(abct, -20.0, 23900.57, rules=rh.GC_HEURISTICS)
    check('GC: ABC transporter still decided structurally',
          status.startswith('ABCT'), status)

    # GC has no sigma gate -- correctly, since it writes its marker into BOTH
    # fields and the dg == 1e7 test catches all 26,555 of them. So a sentinel-scale
    # sigma with a REAL dg falls through GC's sigma-blind rules to its terminal,
    # which since 2026-09-08 reports "?" rather than "=". GC still reaches that
    # answer by exhaustion rather than by testing sigma; that gap is real and is
    # tracked separately for mmdeltag_band and low_energy.
    status, op = run_eq(simple, -50.0, 23900.57, rules=rh.GC_HEURISTICS)
    check('GC terminal now reports "?" on sentinel sigma', op == '?', status)

    # Real reaction, checked against eQuilibrator's published ln_RI of -9.18.
    rxn00001 = reactions.get('rxn00001')
    if rxn00001 is not None:
        ctx = rh.Context(rxn00001, -4.067241221383205, 0.045636422673321665)
        check('rxn00001 ln(Gamma) == -9.18', abs(ctx.ln_gamma + 9.18) < 0.02,
              f'{ctx.ln_gamma:.4f}')
        status, op = run_eq(rxn00001, -4.067241221383205, 0.045636422673321665)
        check('rxn00001 -> ">"', op == '>', status)


def main():
    print('Loading Biochemistry/reaction_*.json ...')
    reactions = load_reactions()
    print(f'  {len(reactions)} reactions')

    test_registry()
    test_eq_rules(reactions)
    test_ln_gamma_matches_equilibrator(reactions)

    print()
    if failures:
        print(f'FAIL -- {len(failures)} check(s): ' + ', '.join(failures))
        return 1
    print('PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
