#!/usr/bin/env python3
"""Classify each compound in _reports/Formula_Conflicts.txt into one of three
categories based on how the source-provided formulas disagree:

  H_only            - sources agree on every heavy-atom count but disagree on H.
                      Typically a protonation-state difference between databases.
  heavy_atom        - sources agree on the element set but disagree on at least
                      one heavy-atom count (skeleton disagreement).
  element_set       - sources disagree on which elements are present (e.g. one
                      source has a halogen the others do not; often indicates a
                      mismapped alias to a different compound).

Reads:  Biochemistry/Structures/_reports/Formula_Conflicts.txt
        (columns: cpd_id, format, stage, formula, charge, source_id, source_db)

Writes: Biochemistry/Structures/_reports/Formula_Conflicts_Classified.txt
        (columns: cpd_id, category, n_variants, formulas)

Prints: summary counts and per-category totals to stdout.

Also emits a second grouping cross-tabulated by "does the compound touch any
priority-scope reaction?" when a priority-compound list is provided via
--priority. This is useful for prioritising remaining curation work against
the v2.0.0 release.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict


def parse_formula(formula):
    """Return {element: count} for a Hill-order formula string.
    Treats 'R' as a first-class element (ModelSEED wildcard convention)."""
    if not formula or formula in ('null', 'NULL'):
        return {}
    out = {}
    for m in re.finditer(r'([A-Z][a-z]?|R|\*)(\d*)', formula):
        elem = m.group(1)
        n = int(m.group(2)) if m.group(2) else 1
        out[elem] = out.get(elem, 0) + n
    return out


def classify(variants):
    """Given a list of (formula, charge, source_id, source_db) tuples for one
    compound, return the conflict category."""
    parsed = [parse_formula(f) for f, _, _, _ in variants]
    heavy_atom_sets = {
        frozenset((el, n) for el, n in p.items() if el != 'H')
        for p in parsed
    }
    if len(heavy_atom_sets) == 1:
        # heavy atoms match; either H differs or the variants are truly
        # identical (which shouldn't happen after dedup but we handle it)
        h_counts = {p.get('H', 0) for p in parsed}
        if len(h_counts) > 1:
            return 'H_only'
        return 'single_variant'
    element_sets = {frozenset(p.keys()) for p in parsed}
    if len(element_sets) > 1:
        return 'element_set'
    return 'heavy_atom'


def load_conflicts(path):
    """Read Formula_Conflicts.txt into {cpd_id: [(formula, charge, source_id,
    source_db), ...]}."""
    by_cpd = defaultdict(list)
    with open(path) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 7:
                continue
            cpd, _fmt, _stage, formula, charge, src_id, src_db = row[:7]
            by_cpd[cpd].append((formula, charge, src_id, src_db))
    return by_cpd


def load_priority(path):
    """One compound ID per line."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument('--conflicts',
                    default='Biochemistry/Structures/_reports/Formula_Conflicts.txt',
                    help='Path to Formula_Conflicts.txt (default: %(default)s)')
    ap.add_argument('--output',
                    default='Biochemistry/Structures/_reports/Formula_Conflicts_Classified.txt',
                    help='Output classified TSV (default: %(default)s)')
    ap.add_argument('--priority', default=None,
                    help='Optional path to a compound-ID-per-line priority list '
                         '(e.g. v7.0 template compounds). Adds a priority-vs-'
                         'not cross-tabulation to the summary output.')
    args = ap.parse_args()

    by_cpd = load_conflicts(args.conflicts)
    priority = load_priority(args.priority)

    # Classify each compound
    classified = {}
    for cpd, variants in by_cpd.items():
        classified[cpd] = classify(variants)

    # Summary counts
    total = len(classified)
    by_cat = defaultdict(list)
    for cpd, cat in classified.items():
        by_cat[cat].append(cpd)

    print(f"Formula conflicts: {total} compounds across {sum(len(v) for v in by_cpd.values())} variant rows")
    print(f"  {'Category':<16} {'Compounds':>10}")
    for cat in ('H_only', 'heavy_atom', 'element_set', 'single_variant'):
        n = len(by_cat.get(cat, []))
        if n:
            print(f"  {cat:<16} {n:>10}")

    if priority is not None:
        print(f"\nPriority-set cross-tabulation ({len(priority)} priority compounds):")
        print(f"  {'Category':<16} {'in priority':>12} {'outside':>10}")
        for cat in ('H_only', 'heavy_atom', 'element_set', 'single_variant'):
            in_pri = sum(1 for cpd in by_cat.get(cat, []) if cpd in priority)
            out_pri = len(by_cat.get(cat, [])) - in_pri
            if in_pri + out_pri:
                print(f"  {cat:<16} {in_pri:>12} {out_pri:>10}")

    # Write classified TSV
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t', lineterminator='\n')
        w.writerow(['cpd_id', 'category', 'n_variants', 'formulas'])
        for cpd in sorted(classified):
            variants = by_cpd[cpd]
            uniq_fc = sorted({(f, c) for f, c, _, _ in variants})
            formulas = ';'.join(f'{f}/{c}' for f, c in uniq_fc)
            w.writerow([cpd, classified[cpd], len(uniq_fc), formulas])

    print(f"\nWrote classified conflicts to {args.output}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
