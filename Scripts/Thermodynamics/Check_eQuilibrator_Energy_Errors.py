#!/usr/bin/env python3
"""Flag eQuilibrator records that are failure reports rather than estimates.

eQuilibrator never raises. When it cannot decompose a reaction it returns a
*number*, and that number reaches us looking like every other number. Two
distinct failures are laundered this way, and both are silent:

1. **The zero return.** ``ComponentContribution.standard_dg`` checks for a
   residual -- a component outside both the reactant and group-contribution
   spans -- and on finding one returns ``Q_(0, "kJ/mol") ± RMSE_inf``, throwing
   the computed mean away. ``standard_dg_prime`` then adds the pH/Mg Legendre
   transform on top **unconditionally**, because it has no residual check of its
   own. So a declined reaction arrives as ``0 + transform``: a plausible,
   nonzero, signed energy whose entire content is how many protons the reaction
   releases at pH 7. Where the transform happens to be nil the zero shows
   through bare, which is the cleanest signature we get.

2. **The sentinel uncertainty.** ``RMSE_inf`` is 1e5 kJ/mol, and multiples of it
   appear when several degrees of freedom are unknown. It is a marker, not an
   error bar. Real sigma across the whole ModelSEED table tops out at 65.35
   kcal/mol against a marker of 23,900.57, so the two regimes are separated by
   an empty gap three orders of magnitude wide -- there is no borderline case to
   get wrong.

A third class is ours, not eQuilibrator's, but it produces the same kind of
confident-looking wrong number and so is checked here too:

3. **The collapsed formula.** ``Retrieve_eQuilibrator_Reactions_Energies.py``
   keys its MetaNetX reaction formula on compound id and drops the compartment.
   Any species appearing on both sides therefore cancels, and eQuilibrator
   scores a *different, smaller* reaction than the one we asked about. The
   published ``ln_reversibility_index`` disagrees with the index recomputed from
   ModelSEED's own stoichiometry, which is how we detect it.

Note this script gives transport reactions no special status. They are checked
by exactly the rules above, and they show up under (3) not because crossing a
membrane is special chemistry but because a formula keyed on compound id loses
the only thing that distinguished the two sides.

Usage:
    ./Check_eQuilibrator_Energy_Errors.py                  # scan the database
    ./Check_eQuilibrator_Energy_Errors.py --table          # scan the raw .tbl too
    ./Check_eQuilibrator_Energy_Errors.py --self-test      # prove the checks fire
    ./Check_eQuilibrator_Energy_Errors.py --tsv out.tsv    # per-reaction detail
    ./Check_eQuilibrator_Energy_Errors.py --strict         # exit 1 if anything is flagged
"""
import argparse
import math
import os
import re
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, 'Libs', 'Python'))

from reversibility_index import ln_reversibility_index   # noqa: E402

# --- Thresholds -----------------------------------------------------------
KJ_PER_KCAL = 4.184

# eQuilibrator's RMSE_inf, the "I declined this reaction" marker.
RMSE_INF_KJ = 1.0e5
RMSE_INF_KCAL = RMSE_INF_KJ / KJ_PER_KCAL              # 23,900.57

# The cut sits in the empty gap between the largest genuine sigma observed
# (65.35 kcal/mol) and the smallest marker value. An order of magnitude below
# the marker and two above the largest real value.
SENTINEL_SIGMA_KCAL = 1.0e4 / KJ_PER_KCAL              # 2,390.06

# Largest sigma we are willing to call an error bar. Between this and the
# sentinel cut is a no-man's-land that should be empty; if anything lands there
# the assumption above has broken and we want to hear about it.
MAX_CREDIBLE_SIGMA_KCAL = 100.0

# ModelSEED's own "no value" marker in the thermodynamics lists.
MSDB_SENTINEL_DG = 10000000

# ln_RI agreement tolerance, matching test_eq_heuristics.py.
LN_RI_TOL_ABS, LN_RI_TOL_REL = 0.10, 0.005

EQ_LABEL = 'eQuilibrator'
EQ_TABLE = os.path.join(REPO_ROOT, 'Biochemistry', 'Thermodynamics',
                        'eQuilibrator', 'MetaNetX_Reaction_Energies.tbl')

PLAIN_MEASUREMENT = re.compile(r'^(-?[\d.]+(?:[eE][-+]?\d+)?)\+/-')

# --- Error codes, most severe first --------------------------------------
CODES = [
    ('SENTINEL_DG',      'dG is ModelSEED\'s 1e7 no-value marker'),
    ('NON_FINITE',       'dG or sigma is NaN/Inf'),
    ('UNDECOMPOSABLE',   'sigma at or past RMSE_inf: eQuilibrator declined the reaction'),
    ('ZERO_RETURN',      'dG is exactly 0.0: the discarded-mean branch, transform and all'),
    ('IMPLAUSIBLE_SIGMA', 'sigma above any credible error bar but below the sentinel cut'),
    ('COLLAPSED_FORMULA', 'published ln_RI disagrees with the index recomputed from our stoichiometry'),
    ('NEGATIVE_SIGMA',   'sigma is negative'),
]
CODE_ORDER = [c for c, _ in CODES]
CODE_HELP = dict(CODES)


def check_record(dg, sigma, stoichiometry=None, published_ln_ri=None):
    """Every error code that applies to one (dG, sigma) record.

    ``dg`` and ``sigma`` are kcal/mol. ``stoichiometry`` and
    ``published_ln_ri`` are optional; supply both to enable the collapsed-formula
    check. Returns a list of codes, empty when the record looks sound.

    This is the whole of the detection logic -- the database scan, the raw-table
    scan and the self-test all call exactly this function, so there is one
    definition of "bad record" and not three.
    """
    codes = []
    if dg is None or sigma is None:
        return ['NON_FINITE']

    dg, sigma = float(dg), float(sigma)

    if not (math.isfinite(dg) and math.isfinite(sigma)):
        return ['NON_FINITE']

    if dg == MSDB_SENTINEL_DG:
        codes.append('SENTINEL_DG')

    if sigma < 0:
        codes.append('NEGATIVE_SIGMA')

    abs_sigma = abs(sigma)
    if abs_sigma >= SENTINEL_SIGMA_KCAL:
        codes.append('UNDECOMPOSABLE')
    elif abs_sigma > MAX_CREDIBLE_SIGMA_KCAL:
        codes.append('IMPLAUSIBLE_SIGMA')

    # The bare zero. Only meaningful when it is not already explained by the
    # sentinel -- but we record it separately either way, because a zero with a
    # *small* sigma is a different and more alarming animal: it means the
    # declined value was propagated without its marker.
    if dg == 0.0 and 'SENTINEL_DG' not in codes:
        codes.append('ZERO_RETURN')

    if stoichiometry is not None and published_ln_ri is not None:
        mine = ln_reversibility_index(stoichiometry, dg)
        if mine is not None:
            tol = max(LN_RI_TOL_ABS, LN_RI_TOL_REL * abs(published_ln_ri))
            if abs(mine - published_ln_ri) > tol:
                codes.append('COLLAPSED_FORMULA')

    return codes


# --- Scanners -------------------------------------------------------------
def scan_database(reactions):
    """Check ``thermodynamics['eQuilibrator']`` for every reaction that has one.

    Reads the per-source dictionary, not the flat ``deltag``/``deltagerr``
    fields: since the additive-thermodynamics refactor nothing rewrites the flat
    fields, so they hold whichever source last promoted to canonical. Checking
    them would be checking Group Contribution's numbers under eQuilibrator's
    name.
    """
    findings = []
    checked = 0
    for rxn_id, entry in reactions.items():
        thermo = entry.get('thermodynamics')
        if not isinstance(thermo, dict) or EQ_LABEL not in thermo:
            continue
        pair = thermo[EQ_LABEL]
        checked += 1
        dg = pair[0] if pair and pair[0] is not None else None
        sigma = pair[1] if len(pair) > 1 and pair[1] is not None else None
        codes = check_record(dg, sigma)
        if codes:
            findings.append((rxn_id, dg, sigma, codes,
                             entry.get('is_transport', 0)))
    return checked, findings


def scan_table(path, reactions=None):
    """Check the raw retrieval output, including the published ln_RI column.

    The database scan cannot see the collapsed-formula class, because column 4
    is discarded on the way in (``_thermo_helpers.parse_two_col_energy_table``
    reads only columns 1 and 2). This is where it surfaces.
    """
    findings = []
    checked = 0
    if not os.path.exists(path):
        return 0, findings
    for line in open(path):
        array = line.rstrip('\n').split('\t')
        if len(array) < 3:
            continue
        rxn_id = array[0]
        checked += 1
        try:
            dg, sigma = float(array[1]), float(array[2])
        except ValueError:
            findings.append((rxn_id, None, None, ['NON_FINITE'], 0))
            continue

        published = None
        stoich = None
        if len(array) > 3 and reactions is not None and rxn_id in reactions:
            match = PLAIN_MEASUREMENT.match(array[3])
            if match:
                published = float(match.group(1))
                stoich = reactions[rxn_id]['stoichiometry']

        codes = check_record(dg, sigma, stoich, published)
        if codes:
            findings.append((rxn_id, dg, sigma, codes,
                             reactions.get(rxn_id, {}).get('is_transport', 0)
                             if reactions else 0))
    return checked, findings


# --- Self-test ------------------------------------------------------------
def self_test():
    """Synthetic records covering each class. The two the docstring promises --
    the 0-return energy and the 1e5 sentinel -- are asserted explicitly."""
    # A -> B, one substrate one product, so sum|nu| = 2.
    simple = [{'compound': 'cpd00100', 'compartment': 0, 'coefficient': -1},
              {'compound': 'cpd00200', 'compartment': 0, 'coefficient': 1}]
    # A[0] -> A[1]: ordinary chemistry, two species, nothing special about it.
    transport = [{'compound': 'cpd00100', 'compartment': 0, 'coefficient': -1},
                 {'compound': 'cpd00100', 'compartment': 1, 'coefficient': 1}]

    cases = [
        # (label, dg, sigma, stoich, published_ln_ri, codes that MUST appear)
        ('the 1e5 kJ/mol sentinel, exactly',
         -20.46, RMSE_INF_KCAL, None, None, {'UNDECOMPOSABLE'}),
        ('a multiple of the sentinel (several unknown DOF)',
         5.0, 33 * RMSE_INF_KCAL, None, None, {'UNDECOMPOSABLE'}),
        ('the sentinel written in kJ/mol by mistake',
         5.0, RMSE_INF_KJ, None, None, {'UNDECOMPOSABLE'}),
        ('the 0-return energy, transform nil, marker intact',
         0.0, RMSE_INF_KCAL, None, None, {'ZERO_RETURN', 'UNDECOMPOSABLE'}),
        ('the 0-return energy with its marker stripped',
         0.0, 0.5, None, None, {'ZERO_RETURN'}),
        ('ModelSEED 1e7 no-value marker',
         MSDB_SENTINEL_DG, 0.0, None, None, {'SENTINEL_DG'}),
        ('NaN dG', float('nan'), 1.0, None, None, {'NON_FINITE'}),
        ('infinite sigma', 1.0, float('inf'), None, None, {'NON_FINITE'}),
        ('missing sigma', 1.0, None, None, None, {'NON_FINITE'}),
        ('negative sigma', 1.0, -2.0, None, None, {'NEGATIVE_SIGMA'}),
        ('sigma in the empty no-man\'s-land',
         1.0, 500.0, None, None, {'IMPLAUSIBLE_SIGMA'}),
        ('collapsed formula: published index double ours',
         -10.0, 1.0, transport, 2 * ln_reversibility_index(transport, -10.0),
         {'COLLAPSED_FORMULA'}),
        # Negative controls: these must produce nothing at all.
        ('a healthy record', -4.067, 0.0456, simple,
         ln_reversibility_index(simple, -4.067), set()),
        ('a healthy record with a large but credible sigma',
         12.0, 65.35, simple, ln_reversibility_index(simple, 12.0), set()),
        ('a healthy transport record, scored as ordinary chemistry',
         -10.0, 1.0, transport, ln_reversibility_index(transport, -10.0), set()),
    ]

    print('Self-test: each synthetic record must produce exactly the '
          'expected codes\n')
    failures = 0
    for label, dg, sigma, stoich, published, expected in cases:
        got = set(check_record(dg, sigma, stoich, published))
        ok = got == expected
        failures += not ok
        shown = ','.join(sorted(got)) or '(clean)'
        want = ','.join(sorted(expected)) or '(clean)'
        print(f"  [{'ok  ' if ok else 'FAIL'}] {label}")
        print(f"           got {shown}" + ('' if ok else f'   expected {want}'))

    # The two headline guarantees, stated as their own assertions so a
    # regression in either is impossible to read past.
    print('\n  Headline guarantees:')
    zero_caught = 'ZERO_RETURN' in check_record(0.0, 0.5)
    sentinel_caught = 'UNDECOMPOSABLE' in check_record(-20.46, RMSE_INF_KCAL)
    for label, ok in (('0-return energy is caught', zero_caught),
                      ('1e5 kJ/mol sentinel sigma is caught', sentinel_caught)):
        failures += not ok
        print(f"  [{'ok  ' if ok else 'FAIL'}] {label}")

    print('\n' + ('PASS' if not failures else f'FAIL ({failures})'))
    return failures


# --- Reporting ------------------------------------------------------------
def summarise(title, checked, findings, out=sys.stdout):
    print(f'\n{title}', file=out)
    print('-' * len(title), file=out)
    print(f'  records checked : {checked}', file=out)
    print(f'  records flagged : {len(findings)}'
          f'  ({100.0 * len(findings) / checked:.1f}%)' if checked else
          '  records flagged : 0', file=out)
    counts = {}
    transport_counts = {}
    for _, _, _, codes, is_transport in findings:
        for code in codes:
            counts[code] = counts.get(code, 0) + 1
            if is_transport == 1:
                transport_counts[code] = transport_counts.get(code, 0) + 1
    if not counts:
        print('  nothing flagged.', file=out)
        return counts
    width = max(len(c) for c in counts)
    print(f'\n  {"code".ljust(width)}  {"count":>7}  {"transport":>9}  meaning',
          file=out)
    for code in CODE_ORDER:
        if code not in counts:
            continue
        print(f'  {code.ljust(width)}  {counts[code]:7d}  '
              f'{transport_counts.get(code, 0):9d}  {CODE_HELP[code]}', file=out)
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--self-test', action='store_true',
                        help='run the synthetic-record checks and exit')
    parser.add_argument('--table', action='store_true',
                        help='also scan the raw MetaNetX_Reaction_Energies.tbl')
    parser.add_argument('--tsv', metavar='PATH',
                        help='write per-reaction findings here')
    parser.add_argument('--strict', action='store_true',
                        help='exit 1 when anything is flagged')
    args = parser.parse_args()

    if args.self_test:
        sys.exit(1 if self_test() else 0)

    from BiochemPy import Reactions
    reactions = Reactions().loadReactions()

    checked, findings = scan_database(reactions)
    summarise(f'Database: thermodynamics["{EQ_LABEL}"]', checked, findings)

    all_findings = [('database',) + f for f in findings]

    if args.table:
        t_checked, t_findings = scan_table(EQ_TABLE, reactions)
        summarise(f'Raw table: {os.path.relpath(EQ_TABLE, REPO_ROOT)}',
                  t_checked, t_findings)
        all_findings += [('table',) + f for f in t_findings]

    if args.tsv:
        with open(args.tsv, 'w') as handle:
            handle.write('source\treaction\tdeltag\tsigma\tis_transport\tcodes\n')
            for source, rxn_id, dg, sigma, codes, is_transport in all_findings:
                handle.write(f'{source}\t{rxn_id}\t{dg}\t{sigma}\t'
                             f'{is_transport}\t{",".join(sorted(codes))}\n')
        print(f'\nWrote {len(all_findings)} findings to {args.tsv}')

    if args.strict and all_findings:
        sys.exit(1)


if __name__ == '__main__':
    main()
