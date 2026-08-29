#!/usr/bin/env python3
"""Rewrite eQuilibrator's failure returns as ModelSEED's "no energy" sentinel.

The two projects both have a way of saying "no value", and they are not
compatible. ModelSEED says it loudly, with a value nothing can mistake for an
energy:

    thermodynamics['<source>'] == [10000000.0, 10000000.0, '?']

eQuilibrator says it quietly, by returning a *plausible number* alongside a
1e5 kJ/mol uncertainty. ``ComponentContribution.standard_dg`` hits a residual --
a component outside both the reactant and group-contribution spans -- and
returns ``Q_(0, "kJ/mol") ± RMSE_inf``, discarding the computed mean;
``standard_dg_prime`` then adds the pH/Mg Legendre transform on top, because it
has no residual check of its own. What lands in the database is
``0 + transform``: signed, nonzero, and carrying no chemistry whatsoever.

Stored side by side, those two conventions mean a consumer that correctly skips
ModelSEED's sentinel will happily consume eQuilibrator's. This script translates
the second into the first.

What counts as a sentinel
-------------------------
**``undecomposable`` (on by default).** ``|sigma| >= 1e4 kJ/mol``. This is the
only unambiguous class, and it is unambiguous in both directions: the residual
branch *always* sets sigma to ``RMSE_inf`` (or a multiple, when several degrees
of freedom are unknown), and the largest genuine sigma in the whole ModelSEED
table is 65.35 kcal/mol against a marker of 23,900.57. The cut sits in an empty
gap three orders of magnitude wide -- measured, not assumed: no record anywhere
in the database has a sigma between 100 and 2,390 kcal/mol.

**``zero_dg`` (OFF by default -- read this before enabling).** ``dG == 0.0``
with a *credible* sigma. It is tempting to treat this as the residual branch
showing through where the transform happened to be nil, and for a handful of
records it is. But the residual branch always stamps ``RMSE_inf`` on the sigma,
so by construction these are **not** that branch. Of the 288 such records, 280
carry ``sigma == 0.0`` exactly and 232 of those are one-substrate/one-product
isomerizations -- L-lysine to D-lysine, D-threo- to D-erythro-isocitrate,
16alpha- to 16beta-hydroxysteroid. eQuilibrator's decomposition is stereo-blind:
both sides decompose to identical groups, so the difference is exactly zero with
exactly zero propagated error. That is a real, if uninformative, output of the
model rather than a failure report, and blanking it discards the (correct)
statement that the two stereoisomers are isoenergetic. Enable this only if you
want that behaviour, and prefer ``--zero-dg-nonzero-sigma`` if what you actually
want is the 8 records where a zero energy is genuinely suspicious.

**``non_finite`` (on by default).** NaN or Inf in either field. Always garbage.

Deliberately *not* included: ``COLLAPSED_FORMULA``. That is a wrong but finite
energy produced by our own retrieval step keying its MetaNetX formula on
compound id and dropping the compartment -- a different defect needing a
different fix (rebuild the formula), not a sentinel. Transport reactions get no
special handling here; they are affected only insofar as that collapse hits them
more often, which this script does not act on.

Usage:
    ./Normalize_eQuilibrator_Sentinels.py                    # dry run, report only
    ./Normalize_eQuilibrator_Sentinels.py --apply            # rewrite the JSON
    ./Normalize_eQuilibrator_Sentinels.py --classes undecomposable non_finite zero_dg
    ./Normalize_eQuilibrator_Sentinels.py --zero-dg-nonzero-sigma
    ./Normalize_eQuilibrator_Sentinels.py --table out.tbl    # normalise the raw table
    ./Normalize_eQuilibrator_Sentinels.py --tsv changes.tsv
    ./Normalize_eQuilibrator_Sentinels.py --self-test

Nothing is written without ``--apply`` (or ``--table``, which only ever writes
the new file you name).
"""
import argparse
import math
import os
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, 'Libs', 'Python'))

# One definition of "bad record", shared with the checker. If the detection
# thresholds ever move, they move for both.
from Check_eQuilibrator_Energy_Errors import (          # noqa: E402
    EQ_LABEL, EQ_TABLE, SENTINEL_SIGMA_KCAL, check_record,
)

# --- ModelSEED's convention ----------------------------------------------
# Matches _thermo_helpers.DEFAULT_DG_DGE, with the direction operator that
# Estimate_Reaction_Reversibility.py assigns to an unusable energy. Verified
# against the 28,689 Group-Contribution records already stored this way: every
# one is exactly [10000000.0, 10000000.0, '?'].
MODELSEED_SENTINEL_DG = 10000000.0
MODELSEED_SENTINEL_OPERATOR = '?'


def modelseed_sentinel(operator=MODELSEED_SENTINEL_OPERATOR):
    """A fresh sentinel triple. Built per call so callers cannot alias one
    shared list into thousands of reaction records."""
    return [MODELSEED_SENTINEL_DG, MODELSEED_SENTINEL_DG, operator]


# --- Which error codes are sentinels -------------------------------------
CLASS_CODES = {
    'undecomposable': ('UNDECOMPOSABLE',),
    'non_finite': ('NON_FINITE',),
    'zero_dg': ('ZERO_RETURN',),
}
DEFAULT_CLASSES = ('undecomposable', 'non_finite')


def to_modelseed_sentinel(dg, sigma, classes=DEFAULT_CLASSES,
                          zero_dg_requires_sigma=False):
    """The MSDB sentinel triple if this record is an eQuilibrator failure
    return, else ``None``.

    This is the whole of the translation. ``zero_dg_requires_sigma`` narrows the
    ``zero_dg`` class to records whose sigma is nonzero, which excludes the
    stereo-blind isomerizations described in the module docstring.
    """
    if dg is not None and float(dg) == MODELSEED_SENTINEL_DG:
        return None                       # already ModelSEED's sentinel

    codes = set(check_record(dg, sigma))
    if not codes:
        return None

    wanted = set()
    for name in classes:
        wanted.update(CLASS_CODES.get(name, ()))

    hit = codes & wanted
    if not hit:
        return None

    if hit == {'ZERO_RETURN'} and zero_dg_requires_sigma:
        try:
            if float(sigma) == 0.0:
                return None
        except (TypeError, ValueError):
            pass

    return modelseed_sentinel()


def _reason(dg, sigma):
    return ','.join(sorted(check_record(dg, sigma))) or '(clean)'


# --- Database ------------------------------------------------------------
def normalize_database(reactions, classes, zero_dg_requires_sigma,
                       label=EQ_LABEL):
    """Translate every sentinel record in ``thermodynamics[label]``.

    Mutates ``reactions`` in memory and returns the change list. Reads and
    writes the per-source dictionary only -- never the flat ``deltag`` /
    ``deltagerr`` fields, which since the additive-thermodynamics refactor carry
    whichever source last promoted to canonical. Blanking those would blank
    another source's number.
    """
    changes = []
    for rxn_id in sorted(reactions):
        entry = reactions[rxn_id]
        thermo = entry.get('thermodynamics')
        if not isinstance(thermo, dict) or label not in thermo:
            continue
        triple = thermo[label]
        dg = triple[0] if triple else None
        sigma = triple[1] if len(triple) > 1 else None

        replacement = to_modelseed_sentinel(
            dg, sigma, classes, zero_dg_requires_sigma)
        if replacement is None:
            continue

        # Preserve nothing from the old triple: the operator was derived from
        # the very number we are discarding.
        changes.append({
            'reaction': rxn_id,
            'old_dg': dg,
            'old_sigma': sigma,
            'old_operator': triple[2] if len(triple) > 2 else '',
            'reason': _reason(dg, sigma),
            'canonical_matches': _canonical_matches(entry, dg),
        })
        thermo[label] = replacement
    return changes


def _canonical_matches(entry, dg):
    """Whether the flat ``deltag`` currently equals the value being discarded.

    Reported, never acted on. If it matches, the canonical energy is itself a
    laundered failure return and needs
    ``Promote_Reaction_Thermodynamics_to_Canonical.py`` re-run -- but deciding
    the canonical value is that script's job, not this one's.
    """
    try:
        return float(entry.get('deltag')) == float(dg)
    except (TypeError, ValueError):
        return False


# --- Raw table -----------------------------------------------------------
def normalize_table(source_path, dest_path, classes, zero_dg_requires_sigma):
    """Write a copy of the retrieval table with sentinels translated.

    Normalising at the table means the sentinel never enters the JSON in the
    first place, which is the better place to fix this if the pipeline is being
    re-run rather than retro-corrected.
    """
    changes = []
    written = 0
    with open(source_path) as src, open(dest_path, 'w') as dest:
        for line in src:
            array = line.rstrip('\n').split('\t')
            if len(array) < 3:
                dest.write(line)
                continue
            try:
                dg, sigma = float(array[1]), float(array[2])
            except ValueError:
                dg = sigma = float('nan')

            replacement = to_modelseed_sentinel(
                dg, sigma, classes, zero_dg_requires_sigma)
            if replacement is not None:
                changes.append({
                    'reaction': array[0], 'old_dg': dg, 'old_sigma': sigma,
                    'old_operator': '', 'reason': _reason(dg, sigma),
                    'canonical_matches': False,
                })
                array[1] = repr(replacement[0])
                array[2] = repr(replacement[1])
                if len(array) > 3:
                    # Column 4 is the published ln_RI, computed from the number
                    # we just discarded.
                    array[3] = MODELSEED_SENTINEL_OPERATOR
            dest.write('\t'.join(array) + '\n')
            written += 1
    return changes, written


# --- Self-test -----------------------------------------------------------
def self_test():
    """Every branch of the translation, including what must be left alone."""
    S = MODELSEED_SENTINEL_DG
    RMSE_INF_KCAL = 1.0e5 / 4.184
    both = ('undecomposable', 'non_finite', 'zero_dg')

    cases = [
        # (label, dg, sigma, classes, requires_sigma, expected triple or None)
        ('the 1e5 sentinel becomes ModelSEED\'s',
         -20.46, RMSE_INF_KCAL, DEFAULT_CLASSES, False, [S, S, '?']),
        ('a multiple of the sentinel too',
         5.0, 33 * RMSE_INF_KCAL, DEFAULT_CLASSES, False, [S, S, '?']),
        ('zero dG carrying the marker',
         0.0, RMSE_INF_KCAL, DEFAULT_CLASSES, False, [S, S, '?']),
        ('NaN', float('nan'), 1.0, DEFAULT_CLASSES, False, [S, S, '?']),
        ('Inf sigma', 1.0, float('inf'), DEFAULT_CLASSES, False, [S, S, '?']),
        ('missing sigma', 1.0, None, DEFAULT_CLASSES, False, [S, S, '?']),

        # zero_dg is off by default and must stay off.
        ('bare zero dG is LEFT ALONE by default',
         0.0, 0.5, DEFAULT_CLASSES, False, None),
        ('bare zero dG is converted when asked',
         0.0, 0.5, both, False, [S, S, '?']),
        ('stereo-blind zero (sigma==0) spared by --zero-dg-nonzero-sigma',
         0.0, 0.0, both, True, None),
        ('suspicious zero (sigma>0) still caught by it',
         0.0, 0.7, both, True, [S, S, '?']),

        # Must never fire.
        ('a healthy record', -4.067, 0.0456, both, False, None),
        ('a large but credible sigma', 12.0, 65.35, both, False, None),
        ('already ModelSEED\'s sentinel is idempotent', S, S, both, False, None),
        ('a genuinely tiny nonzero dG', 1e-9, 0.5, both, False, None),
    ]

    print('Self-test: translation of eQuilibrator sentinels\n')
    failures = 0
    for label, dg, sigma, classes, req, expected in cases:
        got = to_modelseed_sentinel(dg, sigma, classes, req)
        ok = got == expected
        failures += not ok
        print(f"  [{'ok  ' if ok else 'FAIL'}] {label}")
        if not ok:
            print(f'           got {got}, expected {expected}')

    # Structural guarantees.
    print('\n  Structural guarantees:')
    a = to_modelseed_sentinel(0.0, 1e9)
    b = to_modelseed_sentinel(0.0, 1e9)
    checks = [
        ('output matches the stored GC convention exactly',
         a == [10000000.0, 10000000.0, '?']),
        ('each call returns a fresh list (no shared alias)', a is not b),
        ('translating twice is a no-op',
         to_modelseed_sentinel(a[0], a[1]) is None),
        ('the empty gap is real: no sigma between 100 and 2390 is a sentinel',
         to_modelseed_sentinel(5.0, 500.0) is None),
        ('the cut is where the checker says it is',
         to_modelseed_sentinel(5.0, SENTINEL_SIGMA_KCAL) is not None
         and to_modelseed_sentinel(5.0, SENTINEL_SIGMA_KCAL * 0.999) is None),
    ]
    for label, ok in checks:
        failures += not ok
        print(f"  [{'ok  ' if ok else 'FAIL'}] {label}")

    # End-to-end on a synthetic reactions dict.
    print('\n  End-to-end over a reactions dict:')
    fake = {
        'rxn00001': {'deltag': -20.46, 'thermodynamics': {
            EQ_LABEL: [-20.46, RMSE_INF_KCAL, '>'],
            'Group contribution': [-19.0, 3.0, '>']}},
        'rxn00002': {'deltag': -4.067, 'thermodynamics': {
            EQ_LABEL: [-4.067, 0.0456, '>']}},
        'rxn00003': {'deltag': 1.0, 'thermodynamics': {
            'Group contribution': [1.0, 2.0, '=']}},
    }
    changes = normalize_database(fake, DEFAULT_CLASSES, False)
    e2e = [
        ('exactly one record translated', len(changes) == 1),
        ('the right one', changes and changes[0]['reaction'] == 'rxn00001'),
        ('its eQuilibrator triple is now the sentinel',
         fake['rxn00001']['thermodynamics'][EQ_LABEL] == [S, S, '?']),
        ('its Group-Contribution triple is untouched',
         fake['rxn00001']['thermodynamics']['Group contribution'] == [-19.0, 3.0, '>']),
        ('the flat deltag is untouched',
         fake['rxn00001']['deltag'] == -20.46),
        ('canonical collision is reported',
         changes and changes[0]['canonical_matches'] is True),
        ('the healthy record is untouched',
         fake['rxn00002']['thermodynamics'][EQ_LABEL] == [-4.067, 0.0456, '>']),
        ('a reaction with no eQuilibrator record is skipped',
         list(fake['rxn00003']['thermodynamics']) == ['Group contribution']),
        ('a second pass changes nothing',
         normalize_database(fake, DEFAULT_CLASSES, False) == []),
    ]
    for label, ok in e2e:
        failures += not ok
        print(f"  [{'ok  ' if ok else 'FAIL'}] {label}")

    print('\n' + ('PASS' if not failures else f'FAIL ({failures})'))
    return failures


# --- Reporting -----------------------------------------------------------
def report(changes, total_records, applied):
    from collections import Counter
    print(f'\n  eQuilibrator records inspected : {total_records}')
    print(f'  records translated to sentinel : {len(changes)}'
          + (f'  ({100.0 * len(changes) / total_records:.1f}%)'
             if total_records else ''))
    if not changes:
        print('\n  Nothing to do.')
        return
    reasons = Counter(c['reason'] for c in changes)
    width = max(len(r) for r in reasons)
    print(f'\n  {"reason".ljust(width)}  {"count":>7}')
    for reason, count in reasons.most_common():
        print(f'  {reason.ljust(width)}  {count:7d}')

    collisions = sum(1 for c in changes if c['canonical_matches'])
    if collisions:
        print(f'\n  NOTE: {collisions} of these are also the current canonical '
              f'`deltag`.\n  Those flat fields are left alone by design -- re-run '
              f'Promote_Reaction_Thermodynamics_to_Canonical.py\n  to pick a '
              f'replacement from another source.')
    print('\n  ' + ('Written.' if applied else
                    'Dry run -- nothing written. Re-run with --apply.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--apply', action='store_true',
                        help='write the change back to the reactions JSON')
    parser.add_argument('--classes', nargs='+', default=list(DEFAULT_CLASSES),
                        choices=sorted(CLASS_CODES),
                        help='which failure classes count as sentinels '
                             f'(default: {" ".join(DEFAULT_CLASSES)})')
    parser.add_argument('--zero-dg-nonzero-sigma', action='store_true',
                        help='with zero_dg enabled, spare records whose sigma is '
                             'exactly 0 (the stereo-blind isomerizations)')
    parser.add_argument('--source', default=EQ_LABEL,
                        help='thermodynamics dictionary key to normalise')
    parser.add_argument('--table', metavar='DEST',
                        help='normalise the raw retrieval table into DEST '
                             'instead of touching the database')
    parser.add_argument('--tsv', metavar='PATH', help='write the change list here')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()

    if args.self_test:
        sys.exit(1 if self_test() else 0)

    print(f'Sentinel classes : {", ".join(args.classes)}')
    if 'zero_dg' in args.classes:
        print('  zero_dg is enabled -- see the module docstring on the 280 '
              'stereo-blind records')
        if args.zero_dg_nonzero_sigma:
            print('  ...but sigma==0 records are spared')

    if args.table:
        changes, written = normalize_table(
            EQ_TABLE, args.table, args.classes, args.zero_dg_nonzero_sigma)
        print(f'Source           : {os.path.relpath(EQ_TABLE, REPO_ROOT)}')
        print(f'Destination      : {args.table}')
        report(changes, written, applied=True)
    else:
        from BiochemPy import Reactions
        handler = Reactions()
        reactions = handler.loadReactions()
        print(f'Source           : thermodynamics["{args.source}"]  '
              f'(the per-source dictionary, not the flat deltag field)')
        changes = normalize_database(
            reactions, args.classes, args.zero_dg_nonzero_sigma, args.source)
        total = sum(1 for r in reactions.values()
                    if isinstance(r.get('thermodynamics'), dict)
                    and args.source in r['thermodynamics'])
        if args.apply and changes:
            handler.saveReactions(reactions)
        report(changes, total, applied=bool(args.apply and changes))

    if args.tsv:
        with open(args.tsv, 'w') as handle:
            handle.write('reaction\told_deltag\told_sigma\told_operator\t'
                         'reason\tcanonical_deltag_matched\n')
            for c in changes:
                handle.write(f"{c['reaction']}\t{c['old_dg']}\t{c['old_sigma']}\t"
                             f"{c['old_operator']}\t{c['reason']}\t"
                             f"{int(c['canonical_matches'])}\n")
        print(f'\nWrote {len(changes)} rows to {args.tsv}')


if __name__ == '__main__':
    main()
