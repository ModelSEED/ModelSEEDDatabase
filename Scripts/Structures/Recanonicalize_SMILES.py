#!/usr/bin/env python3
"""Mass recanonicalize SMILES strings in the consolidated structure files.

Converts SMILES rows in
  - Biochemistry/Structures/Unique_ModelSEED_Structures.txt  (6-col, header)
  - Biochemistry/Structures/All_ModelSEED_Structures.txt     (8-col, no header)
to RDKit canonical isomeric form via
  Chem.MolToSmiles(MolFromSmiles(s), isomericSmiles=True, canonical=True)

Per-source structure files in Biochemistry/Structures/{KEGG,MetaCyc,ChEBI,Rhea,...}/
are NOT touched -- those preserve the original source strings.

Idempotent: re-running on the output produces no further changes (verified by
the --check-idempotent flag).

Rows where RDKit cannot parse the SMILES are left untouched and counted as
parse failures in the summary. Non-SMILE rows (InChI, InChIKey) are passed
through verbatim.

Usage:
  python Recanonicalize_SMILES.py                  # rewrite both files in place
  python Recanonicalize_SMILES.py --dry-run        # report only, no writes
  python Recanonicalize_SMILES.py --workers 8      # parallelism (default: cpu_count)
  python Recanonicalize_SMILES.py --check-idempotent  # rewrite, then re-run to confirm 0 further changes
"""

import argparse
import csv
import os
import sys
from multiprocessing import Pool, cpu_count

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
STRUCTURES_DIR = os.path.normpath(
    os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry', 'Structures'))

UNIQUE_FILE = os.path.join(STRUCTURES_DIR, 'Unique_ModelSEED_Structures.txt')
ALL_FILE = os.path.join(STRUCTURES_DIR, 'All_ModelSEED_Structures.txt')

# Column layout of each file. Tuple is (smiles_type_value, smiles_col_index_0based,
# type_col_index_0based, has_header).
LAYOUTS = {
    'unique': {
        'path': UNIQUE_FILE,
        'type_col': 1,        # "ID Type Aliases Formula Charge Structure"
        'smiles_col': 5,
        'has_header': True,
    },
    'all': {
        'path': ALL_FILE,
        'type_col': 1,        # "ID Type Stage AliasID Source Formula Charge Structure"
        'smiles_col': 7,
        'has_header': False,
    },
}


def canonical_smiles(smiles):
    """Return (canonical_smiles, status) where status is one of:
       'unchanged', 'updated', 'parse_failure', 'null_or_empty'.

    RDKit's MolToSmiles is not always idempotent for metal-coordinated
    aromatic rings (some Mg/Co/Fe porphyrins oscillate between equivalent
    E/Z bond representations like /C=2) <-> \\C=2)). A few compounds cycle
    through 3+ forms. To guarantee idempotency we run several round-trips
    and pick the lexicographically-smallest form as a deterministic
    tiebreaker. Same molecule in all forms.
    """
    if not smiles or smiles in ('null', 'NULL', '.'):
        return smiles, 'null_or_empty'
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, 'parse_failure'
    # Iterate round-trips until we detect a cycle; then pick the
    # lex-smallest form OF THE CYCLE (excluding any transient lead-in).
    # Transients can appear only on the first run of a freshly-written
    # file, so including them in the min would break idempotency.
    seen = []
    current = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    for _ in range(20):  # ample bound; RDKit cycles in practice are length 1-3
        if current in seen:
            cycle_start = seen.index(current)
            cycle = seen[cycle_start:]
            new = min(cycle)
            break
        seen.append(current)
        nxt_mol = Chem.MolFromSmiles(current)
        if nxt_mol is None:
            new = current
            break
        current = Chem.MolToSmiles(nxt_mol, isomericSmiles=True, canonical=True)
    else:
        # Did not converge in 20 iterations — extreme edge case. Fall back
        # to lex-min of everything seen so far.
        new = min(seen)
    return new, ('updated' if new != smiles else 'unchanged')


def _worker(s):
    return canonical_smiles(s)


def process_file(layout_name, dry_run, workers):
    layout = LAYOUTS[layout_name]
    path = layout['path']
    type_col = layout['type_col']
    smiles_col = layout['smiles_col']
    has_header = layout['has_header']

    with open(path, 'r', newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        rows = list(reader)

    header = rows[0] if has_header else None
    body = rows[1:] if has_header else rows

    # Collect SMILES values to canonicalize, with their row indices
    work = []
    for idx, row in enumerate(body):
        if len(row) <= max(type_col, smiles_col):
            continue
        if row[type_col] != 'SMILE':
            continue
        work.append((idx, row[smiles_col]))

    # Parallel canonicalize
    if workers > 1 and len(work) > 100:
        with Pool(workers) as pool:
            results = pool.map(_worker, [w[1] for w in work], chunksize=512)
    else:
        results = [canonical_smiles(s) for _, s in work]

    counts = {'unchanged': 0, 'updated': 0, 'parse_failure': 0, 'null_or_empty': 0}
    for (idx, _orig), (new_smi, status) in zip(work, results):
        counts[status] += 1
        if status == 'updated':
            body[idx][smiles_col] = new_smi

    if not dry_run:
        # Write back. Use lineterminator='\n' to match Unix convention and avoid
        # introducing CRLF noise that csv default would otherwise emit on Windows.
        with open(path, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
            if has_header:
                writer.writerow(header)
            writer.writerows(body)

    return counts, len(work)


def main():
    parser = argparse.ArgumentParser(
        description="Mass-recanonicalize SMILES in Unique_ and All_ModelSEED_Structures.txt via RDKit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--dry-run', action='store_true',
                        help='Compute changes without writing the files')
    parser.add_argument('--workers', type=int, default=cpu_count(),
                        help='Parallel SMILES processors (default: cpu_count)')
    parser.add_argument('--check-idempotent', action='store_true',
                        help='After rewriting, immediately re-run to confirm 0 further changes')
    parser.add_argument('--only', choices=['unique', 'all', 'both'], default='both',
                        help='Restrict to one file (default: both)')
    args = parser.parse_args()

    targets = ['unique', 'all'] if args.only == 'both' else [args.only]

    for name in targets:
        path = LAYOUTS[name]['path']
        if not os.path.isfile(path):
            print(f"SKIP {name}: not found at {path}", file=sys.stderr)
            continue
        print(f"\n== {name} ({path}) ==")
        counts, total = process_file(name, args.dry_run, args.workers)
        print(f"  SMILE rows seen:    {total}")
        print(f"  unchanged:          {counts['unchanged']}")
        print(f"  updated:            {counts['updated']}")
        print(f"  parse_failure:      {counts['parse_failure']}")
        print(f"  null_or_empty:      {counts['null_or_empty']}")
        if args.dry_run:
            print(f"  (dry-run: file NOT modified)")

    if args.check_idempotent and not args.dry_run:
        print("\n== idempotency check (re-running on the just-written files) ==")
        for name in targets:
            counts, total = process_file(name, dry_run=False, workers=args.workers)
            extra = counts['updated']
            status = 'OK' if extra == 0 else f"FAIL — {extra} additional updates"
            print(f"  {name}: {status}")


if __name__ == '__main__':
    main()
