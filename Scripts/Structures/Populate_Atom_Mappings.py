#!/usr/bin/env python
"""Populate the `atom_mapping` dict on each reaction_*.json from the
outputs of Rebuild_AtomMappings_from_raw.py.

Reads:  Biochemistry/Structures/AtomMappings/all_mapping_no_problem.txt
        Biochemistry/Structures/AtomMappings/rxns_confidence.tsv
Writes: Biochemistry/reaction_*.json  (in place)

Each row of `all_mapping_no_problem.txt` is
`rxnXXXXX cpdAAAAA:E#N=cpdBBBBB:E#M`; rows are grouped by rxn ID and
emitted under the `data` key of a single top-level `atom_mapping` dict:

    "atom_mapping": {
        "data": [
            "cpd00001:O#1=cpd00009:O#2",
            "cpd00012:O#1=cpd00009:O#1",
            ...
        ],
        "confidence": "clean",
        "has_symmetry_groups": false
    }

Field meanings:

- `data` — list of atom-pair strings (`cpdA:E#N=cpdB:E#M`). Endpoints
  may be either a single canonical atom reference or a parenthesised
  set of comma-`;`-joined atom refs when the symmetry-group rewrite
  has been applied (see Build_Atom_Equivalence_Groups.py + future
  extensions to Rebuild_AtomMappings_from_raw.py).

- `confidence` — `clean` (every raw RDT row was already a canonical
  single-pair same-element row) or `salvaged` (at least one raw row was
  a run-on chain / dangling orphan / cross-element pair / malformed,
  and the kept rows are a strict subset). Downstream consumers doing
  mechanism-level tracing (13C flux, exact atom fate) should filter to
  `clean` only.

- `has_symmetry_groups` — true when any row in `data` uses set
  notation for an endpoint. Consumers can use this to enable a
  symmetry-aware renderer (see Solr/README + methods paper) or fall
  back to a strict-single-atom view.

The whole `atom_mapping` dict is added only for reactions present in
the input. Reactions absent from the input have any existing
`atom_mapping` (or the pre-reorganization flat `atom_mapping` /
`atom_mapping_confidence` fields) removed, so re-running after a
refreshed bundle produces a clean state.
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
            old_am = r.get('atom_mapping')
            # Detect the pre-reorganization flat layout so we can migrate
            # cleanly: previously `atom_mapping` was a list of strings and
            # `atom_mapping_confidence` was a sibling top-level string.
            old_flat_conf = r.pop('atom_mapping_confidence', None) if isinstance(old_am, list) else None
            if new_pairs is not None:
                new_conf = confidence.get(rid)
                new_has_sym = any('(' in pair for pair in new_pairs)
                new_dict = {
                    'data': list(new_pairs),
                    'confidence': new_conf or 'clean',
                    'has_symmetry_groups': new_has_sym,
                }
                if old_am != new_dict:
                    r['atom_mapping'] = new_dict
                    if old_am is None:
                        added += 1
                    else:
                        updated += 1
                    file_changed = True
                else:
                    unchanged += 1
                # Migration: strip any lingering flat sibling field.
                if old_flat_conf is not None:
                    file_changed = True
            else:
                if old_am is not None:
                    r.pop('atom_mapping', None)
                    removed += 1
                    file_changed = True
                if old_flat_conf is not None:
                    # already popped above; ensures file write picks up removal
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
