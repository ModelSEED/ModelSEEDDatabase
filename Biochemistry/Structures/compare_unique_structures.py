#!/usr/bin/env python3
"""Compare original and updated Unique_ModelSEED_Structures files.

Produces two CSV reports:
  - Unique_Structures_Changes_Summary.csv: one row per changed compound
    with RDKit-based classification of the change type
  - Unique_Structures_Changes_Statistics.csv: counts per category

Usage:
  python compare_unique_structures.py [original] [updated]

Defaults to Unique_ModelSEED_Structures.txt and
Unique_ModelSEED_Structures_updated.txt in the same directory.
"""

import argparse
import csv
import os
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import FindPotentialStereo, StereoSpecified
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

RDLogger.DisableLog('rdApp.*')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare original and updated Unique structure files.")
    parser.add_argument("original", nargs="?",
                        default=os.path.join(SCRIPT_DIR,
                                             "Unique_ModelSEED_Structures.txt"))
    parser.add_argument("updated", nargs="?",
                        default=os.path.join(SCRIPT_DIR,
                                             "Unique_ModelSEED_Structures_updated.txt"))
    parser.add_argument("--summary", default=os.path.join(
        SCRIPT_DIR, "Unique_Structures_Changes_Summary.csv"))
    parser.add_argument("--statistics", default=os.path.join(
        SCRIPT_DIR, "Unique_Structures_Changes_Statistics.csv"))
    return parser.parse_args()


def load_unique(path):
    data = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        next(reader)
        for row in reader:
            if len(row) < 6:
                continue
            data[(row[0], row[1])] = {
                'aliases': row[2], 'formula': row[3],
                'charge': row[4], 'structure': row[5],
            }
    return data


def stereo_counts(smiles):
    if not smiles:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    info = FindPotentialStereo(mol)
    spec = sum(1 for s in info if s.specified == StereoSpecified.Specified)
    return spec, len(info)


def classify_smile_change(old_smi, new_smi):
    old_mol = Chem.MolFromSmiles(old_smi) if old_smi else None
    new_mol = Chem.MolFromSmiles(new_smi) if new_smi else None
    if not old_mol or not new_mol:
        return 'parse_issue'
    if Chem.MolToSmiles(old_mol) == Chem.MolToSmiles(new_mol):
        return 'canonicalization'
    old_inchi = MolToInchi(old_mol)
    new_inchi = MolToInchi(new_mol)
    old_ik = InchiToInchiKey(old_inchi) if old_inchi else None
    new_ik = InchiToInchiKey(new_inchi) if new_inchi else None
    if not old_ik or not new_ik:
        return 'canonicalization'
    op = old_ik.split('-')
    np_ = new_ik.split('-')
    if old_ik == new_ik:
        return 'smiles_representation_change'
    if op[0] == np_[0] and op[1] == np_[1]:
        return 'protonation_correction'
    if op[0] == np_[0]:
        old_spec, _ = stereo_counts(old_smi)
        new_spec, _ = stereo_counts(new_smi)
        if old_spec is not None and new_spec is not None:
            if new_spec > old_spec:
                return 'stereo_gained'
            if new_spec < old_spec:
                return 'stereo_lost'
        return 'stereo_representation_change'
    return 'connectivity_change'


def main():
    args = parse_args()
    orig = load_unique(args.original)
    updated = load_unique(args.updated)

    all_keys = sorted(set(orig) | set(updated))
    by_cpd = defaultdict(list)
    for key in all_keys:
        cpd_id, typ = key
        o = orig.get(key)
        u = updated.get(key)
        if o is None and u is not None:
            by_cpd[cpd_id].append(('ADDED', typ, None, u))
        elif o is not None and u is None:
            by_cpd[cpd_id].append(('REMOVED', typ, o, None))
        elif o['structure'] != u['structure'] or o['formula'] != u['formula'] \
                or o['charge'] != u['charge'] or o['aliases'] != u['aliases']:
            by_cpd[cpd_id].append(('MODIFIED', typ, o, u))

    summary_rows = []
    for cpd_id in sorted(by_cpd):
        entries = by_cpd[cpd_id]
        change_type = entries[0][0]
        smile_entry = next((e for e in entries if e[1] == 'SMILE'), None)
        if smile_entry is None:
            continue
        _, _, o, u = smile_entry

        if change_type == 'ADDED':
            old_smi, new_smi = '', u['structure']
            old_formula, new_formula = '', u['formula']
            old_charge, new_charge = '', u['charge']
            new_spec, new_pot = stereo_counts(new_smi)
            old_spec, old_pot = '', ''
            if new_spec is None:
                new_spec, new_pot = '', ''
            category = 'added'
        elif change_type == 'REMOVED':
            old_smi, new_smi = o['structure'], ''
            old_formula, new_formula = o['formula'], ''
            old_charge, new_charge = o['charge'], ''
            old_spec, old_pot, new_spec, new_pot = '', '', '', ''
            category = 'removed'
        else:
            old_smi, new_smi = o['structure'], u['structure']
            old_formula, new_formula = o['formula'], u['formula']
            old_charge, new_charge = o['charge'], u['charge']
            old_spec, old_pot = stereo_counts(old_smi)
            new_spec, new_pot = stereo_counts(new_smi)
            if old_spec is None:
                old_spec, old_pot = '', ''
            if new_spec is None:
                new_spec, new_pot = '', ''
            category = classify_smile_change(old_smi, new_smi)

        formula_changed = (old_formula != new_formula
                           and old_formula != '' and new_formula != '')
        charge_changed = (old_charge != new_charge
                          and old_charge != '' and new_charge != '')

        summary_rows.append({
            'cpd_id': cpd_id, 'category': category,
            'old_formula': old_formula, 'new_formula': new_formula,
            'formula_changed': 'yes' if formula_changed else 'no',
            'old_charge': old_charge, 'new_charge': new_charge,
            'charge_changed': 'yes' if charge_changed else 'no',
            'old_stereo_defined': old_spec,
            'new_stereo_defined': new_spec,
            'old_stereo_potential': old_pot,
            'new_stereo_potential': new_pot,
            'old_smiles': old_smi, 'new_smiles': new_smi,
        })

    fieldnames = [
        'cpd_id', 'category',
        'old_formula', 'new_formula', 'formula_changed',
        'old_charge', 'new_charge', 'charge_changed',
        'old_stereo_defined', 'new_stereo_defined',
        'old_stereo_potential', 'new_stereo_potential',
        'old_smiles', 'new_smiles',
    ]
    with open(args.summary, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    cat_counts = Counter(r['category'] for r in summary_rows)
    n_formula = sum(1 for r in summary_rows if r['formula_changed'] == 'yes')
    n_charge = sum(1 for r in summary_rows if r['charge_changed'] == 'yes')

    with open(args.statistics, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['category', 'count'])
        for cat, count in cat_counts.most_common():
            writer.writerow([cat, count])
        writer.writerow(['total_compounds', len(summary_rows)])
        writer.writerow(['formula_changed', n_formula])
        writer.writerow(['charge_changed', n_charge])

    print(f"Summary: {args.summary} ({len(summary_rows)} compounds)")
    print(f"Statistics: {args.statistics}\n")
    for cat, count in cat_counts.most_common():
        print(f"  {cat:<30} {count:>6}")
    print(f"  {'total_compounds':<30} {len(summary_rows):>6}")
    print(f"  {'formula_changed':<30} {n_formula:>6}")
    print(f"  {'charge_changed':<30} {n_charge:>6}")


if __name__ == '__main__':
    main()
