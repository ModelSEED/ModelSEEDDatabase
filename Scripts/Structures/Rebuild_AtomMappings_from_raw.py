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

  1. (withdrawn 2026-09-14 — see "Run-on chains" below)

  2. Two-letter element symbols such as ``Cl``, ``Fe``, ``Mg``, ``Zn``,
     ``Hg`` — the shell regex ``:.#`` matches exactly one character on
     the element position, so any two-letter element blows the pattern.

  3. Dangling half-rows where RDT emits an atom with no partner:

         cpd00011:C#1

  4. All-rows-canonical reactions that are still dropped later by the
     element-pair whitelist stage (which only lists single-letter
     pairs: CC NN OO PP SS BB FF II KK).

Any single occurrence of (2)-(3) in a reaction kills that reaction's
mappings entirely under the shell filter. On the current raw set
(260811 rerun) it accepts 24,267 reactions; this row-level filter keeps
the valid pairs from the rest, reaching ~32,400 reactions.

Run-on chains
=============

This script used to also split run-on chains (``A=B=C=D``) into adjacent
same-element pairs, which took coverage to ~32,900 reactions. That was
wrong and was withdrawn on 2026-09-14.

``run_rdt.sh`` assembles pairs by position, not by identity: it sorts the
per-atom lines by RDT's atom-atom-mapping number and concatenates their
text, discarding the numbers. Every ``from`` line ends ``=`` and every
``to`` line ends ``,``. So whenever RDT leaves an atom unmapped, that
atom's line has no partner and the concatenation glues unrelated atoms
into a chain. A chain is the footprint of atoms RDT *declined to map* --
not a compressed set of true pairs.

Splitting them manufactured 20,292 atom pairs across 4,150 reactions
(1.7% of shipped rows), including 16,512 substrate=substrate and
product=product rows that RDT cannot emit even in principle. Sebastian's
``unite_and_filter_mappings.sh`` rejects all 4,523 chain-carrying
reactions outright; this script kept 4,484 of them.

Recovering the genuine pairs inside a chain requires RDT's AAM numbers,
which all_mapping.txt no longer carries. That needs a fixed run_rdt.sh
and a rerun.

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
    orphan like ``cpd00011:C#1``, nor for the atoms inside a run-on
    chain. Such rows are dropped; other valid pairs in the same
    reaction survive.
  * It does not rewrite RDT's atom indexing — the numeric suffix
    ``#N`` is passed through untouched, so downstream consumers that
    already parse Sebastian's format work unchanged.
  * It does not reject self-mappings (``cpd00001:O#1=cpd00001:O#2``)
    that arrive as genuine single-pair rows, since RDT emits them for
    symmetric intramolecular rearrangements; consumers that only want
    inter-compound edges can filter them at read time. Self-mappings
    fabricated by chain splitting are gone as a side effect of the
    change above.

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
    is *not clean* and yields no pairs at all. Only the reaction's
    genuine single-pair rows survive.

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
    # Chain (3+ pieces) or dangling (0/1 pieces). Both are dropped.
    #
    # A chain is NOT a compressed form of true pairs awaiting unpacking.
    # run_rdt.sh assembles pairs positionally (`sort -n | cut -f3 |
    # tr '\n' ' '`), discarding RDT's atom-atom-mapping numbers, so any
    # atom RDT declined to map leaves its neighbours glued to whatever
    # happened to sort next to them. Splitting `A=B=C=D` into adjacent
    # pairs therefore invents atom-atom relationships that RDT never
    # asserted -- including substrate=substrate rows, which RDT cannot
    # emit at all. Verified 2026-09-14: of 4,523 reactions carrying a
    # chain row, Sebastian's own filter keeps 0; this script used to keep
    # 4,484, contributing 20,292 fabricated pairs (1.7% of shipped rows).
    # See Biochemistry/Structures/AtomMappings/README.md and the forensics
    # writeup for the rxn00010 worked example.
    #
    # Recovering the genuine pairs inside a chain needs the AAM numbers,
    # which are absent from all_mapping.txt -- that requires a fixed
    # run_rdt.sh and an RDT rerun, not a smarter filter here.
    return [], False


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
    n_chain_rows = 0               # dropped, not split -- see module docstring
    chain_rxns = set()
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
            if body.count('=') > 1:
                n_chain_rows += 1
                chain_rxns.add(rxn)
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
    print(f'Run-on chain rows dropped: {n_chain_rows:,} across {len(chain_rxns):,} reactions '
          f'(not split into pairs -- see module docstring)')


if __name__ == '__main__':
    main()
