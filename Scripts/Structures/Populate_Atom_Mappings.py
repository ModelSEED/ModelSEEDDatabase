#!/usr/bin/env python
"""Populate the `atom_mapping` field in reaction_*.json from Sebastian's
UniversalRDT/ModelSEED clean output.

Reads:  Biochemistry/Structures/AtomMappings/all_mapping_no_problem.txt
Writes: Biochemistry/reaction_*.json  (in place)

Each row of the input is `rxnXXXXX cpdAAAAA:E#N=cpdBBBBB:E#M`. Rows are grouped
by rxn ID and the space-delimited-suffix (`cpd:E#N=cpd:E#M`) is emitted as an
entry in a flat list on each reaction.

The field is added only for reactions that appear in the input. Existing
`atom_mapping` fields on reactions absent from the input are removed, so
re-running the script after Sebastian pushes an updated bundle produces a
clean state.
"""
import os
import sys
import json
import glob
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
INPUT = os.path.join(REPO, 'Biochemistry', 'Structures', 'AtomMappings',
                     'all_mapping_no_problem.txt')
RXN_GLOB = os.path.join(REPO, 'Biochemistry', 'reaction_*.json')


def load_mappings(path):
    mappings = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            sp = line.find(' ')
            if sp < 0:
                continue
            rxn_id = line[:sp]
            pair = line[sp + 1:]
            mappings[rxn_id].append(pair)
    return mappings


def main():
    if not os.path.isfile(INPUT):
        sys.exit(f"Input not found: {INPUT}")

    mappings = load_mappings(INPUT)
    sys.stderr.write(f"Loaded {sum(len(v) for v in mappings.values()):,} "
                     f"atom-pair rows for {len(mappings):,} reactions.\n")

    files_written = 0
    added = 0
    updated = 0
    removed = 0
    unchanged = 0
    total_seen = 0

    for path in sorted(glob.glob(RXN_GLOB)):
        with open(path) as fh:
            data = json.load(fh)

        file_changed = False
        for r in data:
            total_seen += 1
            rid = r.get('id')
            new_pairs = mappings.get(rid)
            old_pairs = r.get('atom_mapping')
            if new_pairs is not None:
                if old_pairs is None:
                    r['atom_mapping'] = new_pairs
                    added += 1
                    file_changed = True
                elif old_pairs != new_pairs:
                    r['atom_mapping'] = new_pairs
                    updated += 1
                    file_changed = True
                else:
                    unchanged += 1
            else:
                if old_pairs is not None:
                    r.pop('atom_mapping', None)
                    removed += 1
                    file_changed = True

        if file_changed:
            # sort_keys keeps the on-disk layout stable across runs
            with open(path, 'w') as fh:
                json.dump(data, fh, indent=4, sort_keys=True)
                fh.write('\n')
            files_written += 1

    sys.stderr.write(
        f"\nDone.\n"
        f"  Reactions scanned:  {total_seen:,}\n"
        f"  Added:              {added:,}\n"
        f"  Updated (changed):  {updated:,}\n"
        f"  Removed (stale):    {removed:,}\n"
        f"  Unchanged:          {unchanged:,}\n"
        f"  Files rewritten:    {files_written}\n"
    )


if __name__ == '__main__':
    main()
