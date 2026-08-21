#!/usr/bin/env python
"""Populate the `atom_mapping` and `atom_mapping_confidence` fields in
reaction_*.json from the outputs of Rebuild_AtomMappings_from_raw.py.

Reads:  Biochemistry/Structures/AtomMappings/all_mapping_no_problem.txt
        Biochemistry/Structures/AtomMappings/rxns_confidence.tsv
Writes: Biochemistry/reaction_*.json  (in place)

Each row of `all_mapping_no_problem.txt` is
`rxnXXXXX cpdAAAAA:E#N=cpdBBBBB:E#M`; rows are grouped by rxn ID and the
space-delimited-suffix is emitted as an entry in a flat `atom_mapping`
list on each reaction.

`rxns_confidence.tsv` (produced by Rebuild_AtomMappings_from_raw.py)
tags each mapped reaction as `clean` (every raw RDT row was already
canonical single-pair same-element) or `salvaged` (at least one raw row
was a run-on chain, dangling orphan, cross-element pair, or malformed
and the reaction's kept pairs are a strict subset of the raw output).
The tag is written to the reaction's `atom_mapping_confidence` field so
downstream consumers can filter conservatively where mapping accuracy
matters (e.g., mechanism-level tracing or ¹³C flux analysis).

Both fields are added only for reactions that appear in the input. Any
existing `atom_mapping` / `atom_mapping_confidence` on reactions absent
from the input is removed, so re-running after a refreshed bundle
produces a clean state.
"""
import os
import sys
import json
import glob
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
INPUT = os.path.join(REPO, 'Biochemistry', 'Structures', 'AtomMappings',
                     'all_mapping_no_problem.txt')
CONFIDENCE = os.path.join(REPO, 'Biochemistry', 'Structures', 'AtomMappings',
                          'rxns_confidence.tsv')
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


def load_confidence(path):
    """rxn_id -> "clean" | "salvaged". Empty dict if the file is missing —
    old bundles pre-dating Rebuild_AtomMappings_from_raw.py won't have it,
    and the caller should treat that as "no confidence info available"
    rather than failing hard."""
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as fh:
        header = fh.readline()  # discard header line
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def main():
    if not os.path.isfile(INPUT):
        sys.exit(f"Input not found: {INPUT}")

    mappings = load_mappings(INPUT)
    confidence = load_confidence(CONFIDENCE)
    sys.stderr.write(f"Loaded {sum(len(v) for v in mappings.values()):,} "
                     f"atom-pair rows for {len(mappings):,} reactions.\n")
    if confidence:
        n_c = sum(1 for v in confidence.values() if v == 'clean')
        n_s = sum(1 for v in confidence.values() if v == 'salvaged')
        sys.stderr.write(f"Confidence tags: {n_c:,} clean, {n_s:,} salvaged.\n")
    else:
        sys.stderr.write("(no rxns_confidence.tsv found — atom_mapping_confidence field will not be populated)\n")

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
            new_conf = confidence.get(rid) if new_pairs is not None else None
            old_conf = r.get('atom_mapping_confidence')
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
                if new_conf is not None and new_conf != old_conf:
                    r['atom_mapping_confidence'] = new_conf
                    file_changed = True
                elif new_conf is None and old_conf is not None:
                    r.pop('atom_mapping_confidence', None)
                    file_changed = True
            else:
                if old_pairs is not None:
                    r.pop('atom_mapping', None)
                    removed += 1
                    file_changed = True
                if old_conf is not None:
                    r.pop('atom_mapping_confidence', None)
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
