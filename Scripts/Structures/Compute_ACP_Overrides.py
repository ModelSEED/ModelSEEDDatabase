#!/usr/bin/env python3
"""Compute proposed acps_formula_charge.tsv rows for ACP-named wildcard
compounds not yet in that file, using deltas reverse-engineered from
the existing curated entries.

There are two source-database SMILES conventions in circulation for
ACP compounds; the delta between the raw SMILES formula and the
pantetheine-inclusive override formula is:

  KEGG-style raw SMILES (short, one '*' wildcard, no pantetheine atoms):
      override_formula = raw + C11H21N2O7P     (add phosphopantetheine
                                                  minus the S which is
                                                  already in raw)
      charge = -2                              (two phosphate O deprotonated)

  MetaCyc-style raw SMILES (long, two '*' wildcards, pantetheine already
  included, plus serine backbone atoms and R for either end of the
  protein):
      override_formula = raw - C3H5NO - 1 R    (subtract serine backbone
                                                  and collapse 2R -> 1R)
      charge = -2

The script emits proposed rows to STDOUT (or --output) as a TSV
matching the acps_formula_charge.tsv schema. It NEVER writes to the
override file itself -- the curator reviews the proposal and copies
accepted rows over.

Usage:
  python3 Scripts/Structures/Compute_ACP_Overrides.py
  python3 Scripts/Structures/Compute_ACP_Overrides.py --output proposed.tsv
  python3 Scripts/Structures/Compute_ACP_Overrides.py --name-filter acp
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..'))
OVERRIDES_FILE = os.path.join(
    REPO, 'Biochemistry', 'Curation', 'overrides', 'acps_formula_charge.tsv')
ALL_FILE = os.path.join(
    REPO, 'Biochemistry', 'Structures', 'All_ModelSEED_Structures.txt')
NAMES_FILE = os.path.join(
    REPO, 'Biochemistry', 'Aliases', 'Unique_ModelSEED_Compound_Names.txt')


def load_existing_overrides():
    if not os.path.isfile(OVERRIDES_FILE):
        return set()
    with open(OVERRIDES_FILE) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return {row['ID'].strip() for row in reader if row.get('ID', '').strip()}


def load_names():
    names = {}
    if not os.path.isfile(NAMES_FILE):
        return names
    with open(NAMES_FILE) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0]:
                names.setdefault(parts[0], []).append(parts[1])
    return names


def load_candidate_smiles(name_filter, existing_overrides, names):
    """Return list of dicts for each candidate compound: cpd_id, name,
    representative Charged SMILE, source, source_id.
    Prefers KEGG-style SMILES (short) when both are available since the
    delta computation is more robust for that convention."""
    name_filter_low = name_filter.lower()
    candidates = {}
    with open(ALL_FILE) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 8:
                continue
            cpd, typ, stage, srcid, src, formula, charge, struct = row
            if cpd in existing_overrides:
                continue
            if typ != 'SMILE' or stage != 'Charged':
                continue
            if '*' not in struct:
                continue
            cpd_names = names.get(cpd, [])
            if not any(name_filter_low in n.lower() for n in cpd_names):
                continue
            # Prefer KEGG (shorter SMILES = easier delta); take first-seen otherwise
            existing = candidates.get(cpd)
            if existing is None or (src == 'KEGG' and existing['source'] != 'KEGG'):
                candidates[cpd] = {
                    'cpd_id': cpd, 'name': cpd_names[0] if cpd_names else '',
                    'smiles': struct, 'source': src, 'source_id': srcid,
                }
    return sorted(candidates.values(), key=lambda x: x['cpd_id'])


def compute_override(smiles):
    """Return (formula_dict, charge, style) or (None, None, reason)."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        return None, None, 'rdkit_unavailable'

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, 'smiles_parse_failure'

    atoms = defaultdict(int)
    n_wild = 0
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym == '*':
            n_wild += 1
        else:
            atoms[sym] += 1
    atoms['H'] = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())

    if n_wild == 1:
        atoms['C'] += 11
        atoms['H'] += 21
        atoms['N'] += 2
        atoms['O'] += 7
        atoms['P'] += 1
        atoms['R'] = 1
        return atoms, -2, 'KEGG-style'
    elif n_wild == 2:
        atoms['C'] -= 3
        atoms['H'] -= 5
        atoms['N'] -= 1
        atoms['O'] -= 1
        atoms['R'] = 1
        return atoms, -2, 'MetaCyc-style'
    else:
        return None, None, f'unexpected_wildcard_count_{n_wild}'


def format_formula_hill(atoms):
    """Hill order (C, H, then alphabetical); skip zero counts."""
    parts = []
    for e in ['C', 'H']:
        n = atoms.get(e, 0)
        if n == 0:
            continue
        parts.append(f"{e}{n if n != 1 else ''}")
    for e in sorted(k for k in atoms if k not in ('C', 'H')):
        n = atoms.get(e, 0)
        if n == 0:
            continue
        parts.append(f"{e}{n if n != 1 else ''}")
    return ''.join(parts)


def main():
    parser = argparse.ArgumentParser(
        description="Compute proposed acps_formula_charge.tsv rows",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument('--output', default='-',
                        help='TSV output path (default: stdout)')
    parser.add_argument('--name-filter', default='acp',
                        help='Case-insensitive substring; only compounds whose '
                             'name contains this are considered (default: "acp")')
    args = parser.parse_args()

    existing = load_existing_overrides()
    names = load_names()
    print(f"Existing acps_formula_charge.tsv entries: {len(existing)}",
          file=sys.stderr)

    candidates = load_candidate_smiles(args.name_filter, existing, names)
    print(f"Wildcard candidates matching '{args.name_filter}' "
          f"and not already in overrides: {len(candidates)}", file=sys.stderr)

    out_fh = sys.stdout if args.output == '-' else open(args.output, 'w')
    writer = csv.writer(out_fh, delimiter='\t', lineterminator='\n')
    writer.writerow(['ID', 'name', 'formula', 'charge',
                     'source_style', 'source', 'notes'])

    stats = defaultdict(int)
    for c in candidates:
        atoms, charge, style = compute_override(c['smiles'])
        if atoms is None:
            writer.writerow([c['cpd_id'], c['name'], '', '', '',
                             f"{c['source']}:{c['source_id']}",
                             f'skipped: {style}'])
            stats['skipped'] += 1
        else:
            formula = format_formula_hill(atoms)
            writer.writerow([c['cpd_id'], c['name'], formula, str(charge),
                             style, f"{c['source']}:{c['source_id']}",
                             f'auto-computed via {style} delta'])
            stats[style] += 1

    if out_fh is not sys.stdout:
        out_fh.close()

    print(f"\nSummary:", file=sys.stderr)
    for k, v in sorted(stats.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}", file=sys.stderr)
    print(f"\nProposed rows written to: "
          f"{args.output if args.output != '-' else 'stdout'}", file=sys.stderr)
    print(f"Review the output; append accepted rows to "
          f"{OVERRIDES_FILE}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
