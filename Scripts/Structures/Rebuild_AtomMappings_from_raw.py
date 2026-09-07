#!/usr/bin/env python
"""Rebuild ``all_mapping_no_problem.txt`` from the raw ``all_mapping.txt``
under a row-level filter that recovers per-reaction atom pairs the
shell-based upstream filter throws out wholesale.

Reads:   Biochemistry/Structures/AtomMappings/all_mapping.txt
Writes:  Biochemistry/Structures/AtomMappings/all_mapping_no_problem.txt

Motivation
==========

Sebastian's upstream ``unite_and_filter_mappings.sh`` filters at the
*reaction* level: it flags every reaction that contains at least one
non-canonical row and then discards the entire reaction — every valid
row along with the bad one. The pattern it uses (roughly
``cpd.....:.#[0-9]*=cpd.....:.#[0-9]*``) rejects:

  1. Run-on chains where RDT emits ``A=B=C=D`` when it can't pick a
     single 1-1 mapping for a symmetric group. Common in
     decarboxylations that collapse a carboxyl to CO2:

         cpd00516:C#6=cpd00516:O#1=cpd00516:O#2=cpd00011:O#1

  2. Two-letter element symbols such as ``Cl``, ``Fe``, ``Mg``, ``Zn``,
     ``Hg`` — the shell regex ``:.#`` matches exactly one character on
     the element position, so any two-letter element blows the pattern.

  3. Dangling half-rows where RDT emits an atom with no partner:

         cpd00011:C#1

  4. All-rows-canonical reactions that are still dropped later by the
     element-pair whitelist stage (which only lists single-letter
     pairs: CC NN OO PP SS BB FF II KK).

Any single occurrence of (1)-(3) in a reaction kills that reaction's
mappings entirely. On the current raw set (260811 rerun), the shell
filter accepts 24,267 reactions cleanly; a row-level filter that
recovers what's actually salvageable in the rest brings it to ~32,900
reactions with useful atom-atom mappings.

For the PlantSEED biomass reachability use case this closes six of the
eight remaining gaps (Biotin, Leucine, Lysine, Phosphopantetheine,
Thiamin diphosphate, UDP-Xylose). The two still-unreached — Glucotropaeolin
and Sinalbin — are a PlantSEED curation gap (no producing reaction
exists), not an atom-mapping gap.

Filter rules
============

For each raw row ``rxnXXXXX <body>`` (space-separated), split ``<body>``
at ``=`` into ordered pieces. For every adjacent pair (piece_i, piece_{i+1}):

  * Both pieces must match ``cpd\\d{5}:[A-Za-z]{1,2}#\\d+`` — the
    canonical atom reference (widening the shell filter's single-char
    element slot to 1-2 chars).
  * Both pieces must share the same element symbol. Mixed-element
    pairs like ``C=N``, ``O=S`` — RDT's near-isomorphism false
    matches — get dropped at the row level rather than kill the
    reaction.

Every surviving pair is emitted as its own canonical row
``rxnXXXXX <atomA>=<atomB>``. Per-reaction dedup collapses the many
duplicate pairs RDT emits for symmetric reactants (two O of O2 both
mapped to the same product atom, etc.).

What this deliberately does NOT do:

  * It does not try to *infer* the partner for a dangling
    orphan like ``cpd00011:C#1``. The dangling row is dropped
    silently; other valid pairs in the same reaction survive.
  * It does not rewrite RDT's atom indexing — the numeric suffix
    ``#N`` is passed through untouched, so downstream consumers that
    already parse Sebastian's format work unchanged.
  * It does not reject self-mappings (``cpd00001:O#1=cpd00001:O#2``),
    since RDT emits them for symmetric intramolecular rearrangements
    and they are legitimate mapping data; consumers that only want
    inter-compound edges can filter them at read time.

Idempotency
===========

Running the script twice produces byte-identical output. Rows are sorted
per-reaction to make diffs review-friendly (and this ordering matches
what ``Populate_Atom_Mappings.py`` writes into ``reaction_*.json``).
"""

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import os
import re
import sys
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
AM = os.path.join(REPO, 'Biochemistry', 'Structures', 'AtomMappings')
RAW = os.path.join(AM, 'all_mapping.txt')
OUT = os.path.join(AM, 'all_mapping_no_problem.txt')
RXN_LIST = os.path.join(AM, 'rxns_no_problems.txt')
CONFIDENCE = os.path.join(AM, 'rxns_confidence.tsv')
# Optional — when present, atom-pair rows are rewritten to use set
# notation for atoms in InChI equivalence groups. Produced by
# Build_Atom_Equivalence_Groups.py.
EQUIV_TABLE = os.path.join(AM, 'species_equivalence_groups.tsv')

ATOM = re.compile(r'^cpd\d{5}:([A-Za-z]{1,2})#\d+$')


def load_equivalence_map(path):
    """Return {cpd_id: {atom_label: set_string, ...}, ...}.

    An atom label present as a key indicates the atom is one of an
    equivalence class in that compound; the value is the set-notation
    string covering the whole class (e.g. "(O#1;O#2)"). Consumers can
    look up an atom and, if present, replace the single-atom reference
    with the set. Absent from the map → the atom is unique (no
    equivalence class), pass through unchanged.

    Returns an empty dict if the table doesn't exist; downstream logic
    treats an empty map as "no rewriting". So callers work identically
    whether or not the equivalence work has been run.
    """
    out = defaultdict(dict)
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        fh.readline()  # discard header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            cpd, groups_str = parts[0], parts[1]
            for set_str in groups_str.split(' '):
                inner = set_str.strip('()')
                for atom_label in inner.split(';'):
                    if atom_label:
                        out[cpd][atom_label] = set_str
    return out


def rewrite_atom_ref(atom_ref, equiv_map):
    """Rewrite a single canonical atom ref (`cpdXXXXX:E#N`) via the
    equivalence map. Returns the set-notation form if the atom is in a
    group, otherwise the original string unchanged."""
    cpd, atom = atom_ref.split(':', 1)
    set_str = equiv_map.get(cpd, {}).get(atom)
    if set_str is None:
        return atom_ref
    return f'{cpd}:{set_str}'


def classify_row(body):
    """Return (list_of_pair_strings, is_clean_row) for one raw row body.

    A row is *clean* when it is exactly a canonical single-pair
    ``atomA=atomB`` with both endpoints matching ``cpd\\d{5}:E#N``
    (element 1-2 chars) AND sharing the same element symbol. Anything
    else — run-on chain, dangling orphan, cross-element pair, malformed —
    is *not clean*; the salvage still yields whatever same-element
    pairs can be recovered.

    Sebastian's warning applies here: a reaction that requires salvage on
    any of its rows was one RDT struggled with, and the rows that *look*
    clean may still be subtly wrong. Downstream consumers should treat
    the reaction's whole mapping as lower confidence.
    """
    pieces = body.split('=')
    if len(pieces) == 2:
        a, b = pieces
        ma = ATOM.match(a)
        mb = ATOM.match(b)
        if ma and mb:
            if ma.group(1) == mb.group(1):
                return [f'{a}={b}'], True    # clean canonical pair
            return [], False                 # canonical but cross-element
        return [], False                     # single-pair but malformed
    # Chain (3+ pieces) or dangling (0/1 pieces). Neither is clean;
    # emit whatever same-element adjacent pairs can be salvaged.
    pairs = []
    for a, b in zip(pieces, pieces[1:]):
        ma = ATOM.match(a)
        mb = ATOM.match(b)
        if ma and mb and ma.group(1) == mb.group(1):
            pairs.append(f'{a}={b}')
    return pairs, False


def main():
    if not os.path.isfile(RAW):
        sys.exit(f'Missing raw input: {RAW}')

    equiv_map = load_equivalence_map(EQUIV_TABLE)
    if equiv_map:
        print(f'Loaded equivalence map for {len(equiv_map):,} compounds '
              f'from {os.path.basename(EQUIV_TABLE)} — set-notation rewrite enabled')
    else:
        print(f'No equivalence map at {EQUIV_TABLE} — '
              f'atom-pair rows will not be symmetry-rewritten')

    n_raw_rows = 0
    n_raw_rxns = 0
    per_rxn = defaultdict(set)     # rxn -> {pair_string, ...}
    salvaged = set()               # rxn ids with at least one non-clean raw row
    seen_rxns = set()

    with open(RAW) as fh:
        for line in fh:
            line = line.rstrip('\n')
            sp = line.find(' ')
            if sp <= 0:
                continue
            rxn = line[:sp]
            body = line[sp + 1:]
            n_raw_rows += 1
            if rxn not in seen_rxns:
                seen_rxns.add(rxn)
                n_raw_rxns += 1
            pairs, is_clean_row = classify_row(body)
            if not is_clean_row:
                salvaged.add(rxn)
            for pair in pairs:
                # Rewrite endpoints via the equivalence map. When map is
                # empty this is a no-op; when populated, both endpoints
                # get replaced with set notation if they belong to any
                # equivalence class in their compound.
                a, b = pair.split('=')
                new_a = rewrite_atom_ref(a, equiv_map)
                new_b = rewrite_atom_ref(b, equiv_map)
                per_rxn[rxn].add(f'{new_a}={new_b}')

    n_out_rows = 0
    n_out_rxns = 0
    n_clean = 0
    n_salvaged = 0
    with open(OUT, 'w') as out:
        for rxn in sorted(per_rxn):
            pairs = per_rxn[rxn]
            if not pairs:
                continue
            n_out_rxns += 1
            for pair in sorted(pairs):
                out.write(f'{rxn} {pair}\n')
                n_out_rows += 1

    with open(RXN_LIST, 'w') as out:
        for rxn in sorted(per_rxn):
            if per_rxn[rxn]:
                out.write(f'{rxn}\n')

    with open(CONFIDENCE, 'w') as out:
        out.write('reaction\tatom_mapping_confidence\n')
        for rxn in sorted(per_rxn):
            if not per_rxn[rxn]:
                continue
            level = 'salvaged' if rxn in salvaged else 'clean'
            if level == 'clean': n_clean += 1
            else: n_salvaged += 1
            out.write(f'{rxn}\t{level}\n')

    dropped_rxns = n_raw_rxns - n_out_rxns
    print(f'Read  {n_raw_rows:>10,} rows across {n_raw_rxns:>6,} reactions from {os.path.basename(RAW)}')
    print(f'Wrote {n_out_rows:>10,} rows across {n_out_rxns:>6,} reactions to {os.path.basename(OUT)}')
    print(f'       + wrote {n_out_rxns:>6,} ids to {os.path.basename(RXN_LIST)}')
    print(f'       + wrote {os.path.basename(CONFIDENCE)}: {n_clean:,} clean, {n_salvaged:,} salvaged')
    print(f'Reactions with no surviving row: {dropped_rxns:,} (all their rows were malformed or element-mismatched)')


if __name__ == '__main__':
    main()
