#!/usr/bin/env python3
"""Interactively curate the chosen structure for a ModelSEED compound.

Records the curator's decision in a per-curator TSV at
  Biochemistry/Curation/overrides/structure_picks/<curator>.tsv
which List_ModelSEED_Structures.py reads on its next run to honor
the manual pick instead of falling back to its DB-cascade tiebreaker.

Usage:
  python3 Curate_Structure_Pick.py <cpd_id>
  python3 Curate_Structure_Pick.py <cpd_id> --curator samseaver
  python3 Curate_Structure_Pick.py <cpd_id> --rationale "matches Lit ref ..."

Prompts for:
  - Which source's structure to pick (or custom paste)
  - Rationale (one line of free text)
"""

import argparse
import csv
import datetime
import getpass
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..'))
ALL_FILE = os.path.join(REPO_ROOT, 'Biochemistry', 'Structures',
                        'All_ModelSEED_Structures.txt')
UNIQUE_FILE = os.path.join(REPO_ROOT, 'Biochemistry', 'Structures',
                           'Unique_ModelSEED_Structures.txt')
PICK_REASONS_FILE = os.path.join(REPO_ROOT, 'Biochemistry', 'Structures',
                                 'Pick_Reasons.txt')
NAMES_FILE = os.path.join(REPO_ROOT, 'Biochemistry', 'Aliases',
                          'Unique_ModelSEED_Compound_Names.txt')
OVERRIDES_DIR = os.path.join(REPO_ROOT, 'Biochemistry', 'Curation',
                             'overrides', 'structure_picks')

HEADER = ['cpd_id', 'format', 'structure', 'source_db', 'source_id',
          'date', 'rationale']


def load_all_for_compound(cpd_id):
    """Return list of dicts: candidate source rows for the compound,
    filtered to InChI/Charged rows (the picker's preference order)."""
    candidates = []
    with open(ALL_FILE) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 8 or row[0] != cpd_id:
                continue
            # row: cpd_id, type, stage, source_id, source_db, formula, charge, structure
            candidates.append({
                'type': row[1], 'stage': row[2], 'source_id': row[3],
                'source_db': row[4], 'formula': row[5], 'charge': row[6],
                'structure': row[7],
            })
    return candidates


def load_unique_for_compound(cpd_id):
    """Return dict {type: structure} from Unique file."""
    result = {}
    with open(UNIQUE_FILE) as fh:
        next(fh)
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 6 or row[0] != cpd_id:
                continue
            result[row[1]] = {'structure': row[5], 'formula': row[3],
                              'charge': row[4]}
    return result


def load_current_pick_reason(cpd_id):
    if not os.path.isfile(PICK_REASONS_FILE):
        return None
    with open(PICK_REASONS_FILE) as fh:
        next(fh)
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 4 or row[0] != cpd_id:
                continue
            return row[3]
    return None


def load_names(cpd_id, max_n=3):
    if not os.path.isfile(NAMES_FILE):
        return []
    names = []
    with open(NAMES_FILE) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0] == cpd_id:
                names.append(parts[1])
                if len(names) >= max_n:
                    break
    return names


def append_pick(curator_file, row):
    """Append a row to the curator's TSV. Creates the file with header
    if it doesn't exist. Returns True if appended, False if user
    declined to overwrite an existing entry."""
    os.makedirs(os.path.dirname(curator_file), exist_ok=True)
    existing_rows = []
    if os.path.isfile(curator_file):
        with open(curator_file) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            existing_rows = list(reader)
        for er in existing_rows:
            if er.get('cpd_id') == row['cpd_id']:
                print(f"\nNote: {row['cpd_id']} already curated by you on "
                      f"{er.get('date')}.")
                print(f"  Existing rationale: {er.get('rationale')}")
                resp = input("Overwrite with new pick? [y/N]: ").strip().lower()
                if resp != 'y':
                    return False
                existing_rows = [r for r in existing_rows
                                 if r.get('cpd_id') != row['cpd_id']]
                break
    existing_rows.append(row)
    with open(curator_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=HEADER, delimiter='\t',
                                lineterminator='\n')
        writer.writeheader()
        for r in existing_rows:
            writer.writerow({k: r.get(k, '') for k in HEADER})
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Interactively curate a ModelSEED compound's structure pick")
    parser.add_argument('cpd_id', help='ModelSEED compound ID (e.g. cpd02168)')
    parser.add_argument('--curator', default=None,
                        help='Curator handle (defaults to $USER)')
    parser.add_argument('--rationale', default=None,
                        help='Rationale (prompted interactively if omitted)')
    args = parser.parse_args()

    cpd_id = args.cpd_id
    curator = args.curator or getpass.getuser()
    today = datetime.date.today().isoformat()

    # Show current state
    print(f"\nCompound: {cpd_id}")
    names = load_names(cpd_id)
    if names:
        print(f"Names: {' | '.join(names)}")
    current = load_unique_for_compound(cpd_id)
    reason = load_current_pick_reason(cpd_id)
    if current:
        print(f"Currently in Unique:")
        for typ in ('InChI', 'InChIKey', 'SMILE'):
            if typ in current:
                s = current[typ]['structure']
                disp = s if len(s) <= 90 else s[:87] + '...'
                print(f"  {typ:8}  {disp}")
        if reason:
            print(f"  Pick reason: {reason}")
    else:
        print(f"  (not currently in Unique file)")

    # Show candidate sources, grouped by (source_db, source_id) pair
    candidates = load_all_for_compound(cpd_id)
    if not candidates:
        print(f"\nNo source rows found for {cpd_id} in All file.")
        return 1

    # We curate at the InChI/Charged level (picker's preferred type/stage);
    # fall back to SMILE/Charged if no InChI Charged row exists.
    preferred = [c for c in candidates if c['type'] == 'InChI'
                 and c['stage'] == 'Charged']
    if not preferred:
        preferred = [c for c in candidates if c['type'] == 'SMILE'
                     and c['stage'] == 'Charged']
    if not preferred:
        print(f"\nNo InChI or SMILE Charged-stage rows for {cpd_id}.")
        return 1

    print(f"\nCandidate {preferred[0]['type']}/Charged structures:")
    for i, c in enumerate(preferred, 1):
        s = c['structure']
        disp = s if len(s) <= 70 else s[:67] + '...'
        print(f"  [{i}] {c['source_db']:8} {c['source_id']:12} {disp}")
    print(f"  [c] custom (paste your own InChI or SMILES)")
    print(f"  [q] quit without recording")

    choice = input("\nChoice: ").strip().lower()
    if choice == 'q':
        print("Aborted.")
        return 0

    if choice == 'c':
        fmt = input("Format [InChI/SMILE]: ").strip()
        if fmt not in ('InChI', 'SMILE'):
            print(f"Invalid format: {fmt}")
            return 1
        struct = input(f"Paste {fmt}: ").strip()
        if not struct:
            print("Empty structure; aborted.")
            return 1
        source_db = 'custom'
        source_id = 'custom'
    else:
        try:
            idx = int(choice) - 1
            chosen = preferred[idx]
        except (ValueError, IndexError):
            print(f"Invalid choice: {choice}")
            return 1
        fmt = chosen['type']
        struct = chosen['structure']
        source_db = chosen['source_db']
        source_id = chosen['source_id']

    rationale = args.rationale
    if rationale is None:
        rationale = input("Rationale (one line): ").strip()
    if not rationale:
        print("Rationale required; aborted.")
        return 1

    row = {
        'cpd_id': cpd_id,
        'format': fmt,
        'structure': struct,
        'source_db': source_db,
        'source_id': source_id,
        'date': today,
        'rationale': rationale,
    }
    curator_file = os.path.join(OVERRIDES_DIR, f'{curator}.tsv')
    if append_pick(curator_file, row):
        print(f"\nWrote pick to {curator_file}:")
        print(f"  {cpd_id}  {fmt}  {source_db}:{source_id}  {rationale[:60]}")
        print(f"\nNext: re-run Scripts/Structures/List_ModelSEED_Structures.py")
        print(f"      to regenerate Unique_ModelSEED_Structures.txt with this pick honored.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
