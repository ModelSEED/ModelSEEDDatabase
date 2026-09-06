#!/usr/bin/env python3
"""Run two or more reversibility rule sets over the *same* energies and diff them.

The question this answers is narrow on purpose: **holding the thermodynamic
input fixed, how much of a direction call is the rule set?** Comparing the
shipped GC directions against the shipped eQuilibrator directions cannot answer
it, because those two runs differ in their energies as well as their rules.
Here every rule set is handed the identical ``(dG, sigma)`` pair, so any
disagreement is attributable to the rules alone.

Energies come from ``thermodynamics[<source>]`` -- the per-source dictionary.
Not the flat ``deltag``/``deltagerr`` fields: since the additive-thermodynamics
refactor nothing rewrites those, so they carry whichever source last promoted to
canonical. Reading them would silently score Group Contribution's numbers under
eQuilibrator's name, which is the defect that motivated this whole comparison.

Transport reactions are ordinary reactions here. The GC and EQ sets contain
structural shortcuts that fire on them before any energy is read, and those
shortcuts are reported as their own rows -- one of the things worth seeing is
how often the shortcut, not the thermodynamics, is what produced the answer.
The ``RI`` set has no shortcuts at all and serves as the control.

Usage:
    ./Compare_Reversibility_Heuristics.py                       # eQuilibrator: GC vs EQ
    ./Compare_Reversibility_Heuristics.py --sets GC EQ EQ2 RI
    ./Compare_Reversibility_Heuristics.py --source dGPredictor --sets GC DGP
    ./Compare_Reversibility_Heuristics.py --tsv out.tsv
"""
import argparse
import os
import sys
from collections import Counter, OrderedDict

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, 'Libs', 'Python'))

import reversibility_heuristics as rh     # noqa: E402

OPERATORS = ('>', '<', '=', '?')


def rule_family(status):
    """Collapse a status string down to the rule that produced it, so the
    per-rule tallies stay readable."""
    if status is None:
        return '(none)'
    for prefix in ('MdeltaG(Max)', 'MdeltaG(Min)', 'mMdeltaG', 'ABCT',
                   'ATPS', 'Incomplete', 'lowE'):
        if status.startswith(prefix):
            return prefix
    if status.startswith(('EQ:', 'RI:')):
        # "EQ:lnGamma: -9.18+/-0.05" -> "EQ:lnGamma"; "EQ:transport-uncorrected"
        # has no payload and is already its own family.
        return ':'.join(status.split(':', 2)[:2])
    return status.split(':', 1)[0] if ':' in status else status


def compare(reactions, source, set_names):
    """Score every reaction carrying ``source`` energies under each rule set."""
    rows = []
    for rxn_id in sorted(reactions):
        entry = reactions[rxn_id]
        pair = rh._thermo_pair(entry, source)
        if pair is None:
            continue
        dg, sigma = pair
        result = OrderedDict()
        for name in set_names:
            status, operator, _ = rh.run_reversibility(
                entry, rh.explicit_energy(dg, sigma), rh.get_heuristics(name))
            result[name] = (status, operator)
        rows.append((rxn_id, dg, sigma, entry.get('is_transport', 0), result))
    return rows


def print_matrix(rows, left, right):
    """Confusion matrix of two rule sets' operators."""
    counts = Counter((r[4][left][1], r[4][right][1]) for r in rows)
    total = len(rows)
    agree = sum(v for (a, b), v in counts.items() if a == b)

    print(f'\n{left} vs {right} -- operator confusion matrix')
    print(f"  {'':>10s}" + ''.join(f'{op:>9s}' for op in OPERATORS)
          + f"{'total':>9s}    <- {right}")
    for a in OPERATORS:
        row = [counts.get((a, b), 0) for b in OPERATORS]
        if not sum(row):
            continue
        print(f'  {a:>10s}' + ''.join(f'{v:9d}' for v in row)
              + f'{sum(row):9d}')
    col_totals = [sum(counts.get((a, b), 0) for a in OPERATORS)
                  for b in OPERATORS]
    print(f"  {'total':>10s}" + ''.join(f'{v:9d}' for v in col_totals)
          + f'{total:9d}')
    print(f'  ^-- {left}')
    print(f'\n  agree {agree}/{total} ({100.0 * agree / total:.1f}%), '
          f'differ {total - agree} ({100.0 * (total - agree) / total:.1f}%)')

    # Where they differ, which way did each go?
    diffs = Counter((a, b) for (a, b), v in counts.items() if a != b
                    for _ in range(v))
    if diffs:
        print('\n  biggest disagreements:')
        for (a, b), v in diffs.most_common(6):
            print(f'    {left} "{a}" -> {right} "{b}" : {v}')


def print_rule_usage(rows, name):
    """Which rule actually fired, and what each one decided."""
    fired = Counter()
    by_op = {}
    for _, _, _, _, result in rows:
        family = rule_family(result[name][0])
        fired[family] += 1
        by_op.setdefault(family, Counter())[result[name][1]] += 1

    print(f'\n{name} -- which rule fired')
    width = max(len(f) for f in fired)
    print(f"  {'rule'.ljust(width)}  {'count':>7s}  " +
          ''.join(f'{op:>7s}' for op in OPERATORS))
    for family, count in fired.most_common():
        ops = by_op[family]
        print(f'  {family.ljust(width)}  {count:7d}  ' +
              ''.join(f'{ops.get(op, 0):7d}' for op in OPERATORS))


def print_transport_split(rows, left, right):
    """The same agreement figure, split by whether the reaction crosses a
    membrane. Transport gets no special scoring here -- this is a diagnostic,
    to show where the two rule sets' structural shortcuts are doing the work."""
    print(f'\n{left} vs {right} -- agreement split by transport')
    for label, want in (('non-transport', 0), ('transport', 1)):
        subset = [r for r in rows if r[3] == want]
        if not subset:
            continue
        agree = sum(1 for r in subset if r[4][left][1] == r[4][right][1])
        print(f'  {label:14s} {agree:6d}/{len(subset):<6d} '
              f'({100.0 * agree / len(subset):5.1f}% agree)')


ALL_SOURCES = ['Group contribution', 'eQuilibrator',
               'dGPredictor', 'dGPredictor-ModelSEED']

SHORT_NAME = {'Group contribution': 'GroupContrib', 'eQuilibrator': 'eQuilibrator',
              'dGPredictor': 'dGPredictor', 'dGPredictor-ModelSEED': 'dGP-ModelSEED'}


def score_all_sources(reactions, set_for_source):
    """``{source: {rxn_id: operator}}``, each source scored with its own rule set.

    ``set_for_source`` maps a ``thermodynamics`` key to a rule-set name, so the
    matrix can mix rule sets. That is the realistic configuration: Group
    contribution keeps the historical cascade it was designed around and is the
    byte-compare anchor, while the prediction sources move to the index. Passing
    the same name for every source gives the uniform-rules matrix instead.

    Every source is scored from its *own* ``thermodynamics[source]`` energy, so
    the direction reflects that prediction and not whichever one last promoted
    to canonical.
    """
    scored = {src: {} for src in set_for_source}
    heuristics = {src: rh.get_heuristics(name)
                  for src, name in set_for_source.items()}
    for rxn_id, entry in reactions.items():
        for src in set_for_source:
            pair = rh._thermo_pair(entry, src)
            if pair is None:
                continue
            _, operator, _ = rh.run_reversibility(
                entry, rh.explicit_energy(*pair), heuristics[src])
            scored[src][rxn_id] = operator
    return scored


def similarity_matrix(scored, sources):
    """Pairwise agreement, computed only over reactions both sources predict.

    Returns ``(agreement, overlap)``: percent identical operators, and the size
    of the intersection the percent is over. The second matters as much as the
    first -- two sources can agree 95% on the 300 reactions they share and say
    nothing about the other 30,000.
    """
    agreement, overlap = {}, {}
    for a in sources:
        for b in sources:
            shared = scored[a].keys() & scored[b].keys()
            overlap[(a, b)] = len(shared)
            same = sum(1 for r in shared if scored[a][r] == scored[b][r])
            agreement[(a, b)] = 100.0 * same / len(shared) if shared else float('nan')
    return agreement, overlap


def print_similarity(scored, sources, set_for_source, title):
    agreement, overlap = similarity_matrix(scored, sources)
    labels = [SHORT_NAME.get(s, s) for s in sources]
    width = max(len(l) for l in labels)

    print(f'\n{title} -- % identical direction')
    print('  rule set per source: '
          + ', '.join(f'{SHORT_NAME.get(s, s)}={set_for_source[s]}'
                      for s in sources))
    print(f"  {'':{width}s}" + ''.join(f'{l:>15s}' for l in labels))
    for a, la in zip(sources, labels):
        cells = ''.join(
            f'{agreement[(a, b)]:14.1f}%' if a != b else f"{'--':>15s}"
            for b in sources)
        print(f'  {la:{width}s}{cells}')

    print(f'\n  ...over this many shared reactions')
    print(f"  {'':{width}s}" + ''.join(f'{l:>15s}' for l in labels))
    for a, la in zip(sources, labels):
        cells = ''.join(f'{overlap[(a, b)]:15d}' for b in sources)
        print(f'  {la:{width}s}{cells}')

    print(f'\n  directions produced per source (own coverage)')
    for src, label in zip(sources, labels):
        ops = Counter(scored[src].values())
        total = len(scored[src])
        hard = ops['>'] + ops['<']
        print(f'  {label:{width}s} n={total:6d}  '
              f'>{ops[">"]:6d}  <{ops["<"]:6d}  ={ops["="]:6d}  ?{ops["?"]:5d}'
              f'   hard {100.0 * hard / total:.1f}%')
    return agreement, overlap


def print_agreement_decomposition(scored, sources, title):
    """Split each pair's agreement into *why* they agree.

    A high agreement number means little on its own: two sources that both
    default to "reversible" agree without either having said anything. Splitting
    the agreement into "both abstained" and "both committed to the same
    direction" -- and reporting outright directional conflicts alongside --
    separates consensus from mutual silence.

    ``excl "?"`` repeats the agreement over reactions where *both* sources made
    a call, so the undecomposable records eQuilibrator correctly refuses do not
    read as disagreement.
    """
    print(f'\n{title} -- what the agreement is made of '
          f'(% of shared reactions)')
    header = (f"  {'pair':30s}{'agree':>8s}{'both =':>9s}{'both ->':>9s}"
              f"{'conflict':>10s}{'excl ?':>9s}")
    print(header)
    for i, a in enumerate(sources):
        for b in sources[i + 1:]:
            shared = scored[a].keys() & scored[b].keys()
            n = len(shared)
            if not n:
                continue
            agree = sum(1 for r in shared if scored[a][r] == scored[b][r])
            both_eq = sum(1 for r in shared
                          if scored[a][r] == '=' and scored[b][r] == '=')
            both_dir = sum(1 for r in shared
                           if scored[a][r] in '<>' and scored[b][r] in '<>'
                           and scored[a][r] == scored[b][r])
            conflict = sum(1 for r in shared
                           if scored[a][r] in '<>' and scored[b][r] in '<>'
                           and scored[a][r] != scored[b][r])
            called = [r for r in shared
                      if scored[a][r] != '?' and scored[b][r] != '?']
            excl = (100.0 * sum(1 for r in called if scored[a][r] == scored[b][r])
                    / len(called)) if called else float('nan')
            label = f'{SHORT_NAME.get(a, a)} / {SHORT_NAME.get(b, b)}'
            print(f'  {label:30s}{100.0 * agree / n:7.1f}%{100.0 * both_eq / n:8.1f}%'
                  f'{100.0 * both_dir / n:8.1f}%{100.0 * conflict / n:9.1f}%'
                  f'{excl:8.1f}%')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--source', default='eQuilibrator',
                        help='thermodynamics dictionary key supplying the energies')
    parser.add_argument('--sets', nargs='+', default=['GC', 'EQ'],
                        help='rule sets to run (%s)' % ', '.join(sorted(rh.HEURISTIC_SETS)))
    parser.add_argument('--tsv', metavar='PATH', help='per-reaction detail')
    parser.add_argument('--matrix', action='store_true',
                        help='cross-source similarity matrix: hold the rule set '
                             'fixed and vary which prediction supplies the energy')
    parser.add_argument('--sources', nargs='+', default=ALL_SOURCES,
                        help='sources for --matrix (default: all four)')
    parser.add_argument('--index-sources', nargs='+',
                        default=['eQuilibrator', 'dGPredictor',
                                 'dGPredictor-ModelSEED'],
                        help='in --matrix, the sources that take the second and '
                             'later --sets rule sets. Everything else keeps the '
                             'first one, so Group contribution stays on the '
                             'cascade it was designed around.')
    args = parser.parse_args()

    for name in args.sets:
        if name not in rh.HEURISTIC_SETS:
            sys.exit(f'ERROR: unknown rule set {name!r}; choose from '
                     + ', '.join(sorted(rh.HEURISTIC_SETS)))
    if not args.matrix and len(args.sets) < 2:
        sys.exit('ERROR: need at least two rule sets to compare')

    from BiochemPy import Reactions
    reactions = Reactions().loadReactions()

    if args.matrix:
        baseline = args.sets[0]
        # First pass: the baseline rule set everywhere -- what ModelSEED does now.
        plans = [({src: baseline for src in args.sources},
                  f'{baseline} rules on every source')]
        # Then one pass per alternate set, applied only to --index-sources.
        for name in args.sets[1:]:
            plans.append((
                {src: (name if src in args.index_sources else baseline)
                 for src in args.sources},
                f'{name} on {", ".join(SHORT_NAME.get(s, s) for s in args.sources if s in args.index_sources)}'
                f'; {baseline} elsewhere'))
        for set_for_source, title in plans:
            scored = score_all_sources(reactions, set_for_source)
            print_similarity(scored, args.sources, set_for_source, title)
            print_agreement_decomposition(scored, args.sources, title)
        return

    rows = compare(reactions, args.source, args.sets)
    print(f'Source   : thermodynamics["{args.source}"]  '
          f'(the per-source dictionary, not the flat deltag field)')
    print(f'Rule sets: {", ".join(args.sets)}')
    print(f'Reactions: {len(rows)}')

    for name in args.sets:
        print_rule_usage(rows, name)

    left = args.sets[0]
    for right in args.sets[1:]:
        print_matrix(rows, left, right)
        print_transport_split(rows, left, right)

    if args.tsv:
        with open(args.tsv, 'w') as handle:
            header = ['reaction', 'deltag', 'sigma', 'is_transport']
            for name in args.sets:
                header += [f'{name}_operator', f'{name}_rule', f'{name}_status']
            handle.write('\t'.join(header) + '\n')
            for rxn_id, dg, sigma, is_transport, result in rows:
                fields = [rxn_id, f'{dg:.6f}', f'{sigma:.6f}', str(is_transport)]
                for name in args.sets:
                    status, operator = result[name]
                    fields += [operator, rule_family(status), status or '']
                handle.write('\t'.join(fields) + '\n')
        print(f'\nWrote {len(rows)} rows to {args.tsv}')


if __name__ == '__main__':
    main()
