#!/usr/bin/env python3
"""Compile the biochemistry JSON into the legacy Solr payload shape
that matches the pre-2026 (master-branch) schema.

Reads:   ../../Biochemistry/compound_*.json
         ../../Biochemistry/reaction_*.json
Writes:  ./solr_compounds_legacy.json
         ./solr_reactions_legacy.json

Purpose
=======

The 2026 SOLR upgrade introduces a substantially expanded schema —
nested per-source thermodynamics, atom_mapping as a structured dict,
stoichiometry as child docs, denormalized flat filters. That schema is
what powers the *staging* UI while it evolves. The *production* UI
still expects the old flat schema (single `deltag` / `deltagerr`
scalars per compound and reaction, `stoichiometry` as a single string,
no atom_mapping / thermodynamics / pkas dicts).

Loading data compiled for the new schema into production cores would
break the current production UI. So this script produces a compatibility
payload matching the pre-upgrade layout, for posting into the
`compounds_legacy` / `reactions_legacy` cores in the same Solr instance.

Contract vs the production (master-branch) schema
=================================================

Compounds — keep flat scalar fields exactly as master's schema declares
them; drop `thermodynamics` / `pkas` dicts (unindexed).

Reactions — keep flat scalar fields; drop `thermodynamics` /
`atom_mapping` dicts; flatten `stoichiometry` (list of dicts) into a
single semicolon-joined string of `coef:cpd:compartment:0:"name"`
tuples — the format master's compounds.json / reactions.json ships.

See the companion Compile_Biochemistry_for_SOLR.py for the *staging*
(new nested) payload. Both scripts read the same source JSONs; only
the output shape differs.
"""
import glob
import json
import os

BIOCHEM_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            '..', '..', 'Biochemistry'))
OUT_COMPOUNDS = os.path.join(os.path.dirname(__file__),
                             'solr_compounds_legacy.json')
OUT_REACTIONS = os.path.join(os.path.dirname(__file__),
                             'solr_reactions_legacy.json')

# Fields that were added to the JSON schema after the master release
# but aren't indexed by the legacy Solr schema. Drop them from the
# compile output so we don't ship data the schema won't consume.
DROP_FROM_COMPOUND = ('thermodynamics', 'pkas', 'class')
DROP_FROM_REACTION = ('thermodynamics', 'atom_mapping')


def flatten_stoichiometry(entries):
    """Convert the reaction's list-of-dicts stoichiometry into the flat
    string master's Solr schema expects:
        -1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"
    Format is `coefficient:compound:compartment:0:"name"`, semicolon-joined.
    Returns '' when the reaction has no stoichiometry entries."""
    parts = []
    for e in entries or ():
        if not isinstance(e, dict):
            continue
        parts.append(':'.join([
            str(e.get('coefficient', 0)),
            str(e.get('compound', '')),
            str(e.get('compartment', 0)),
            '0',
            '"{}"'.format((e.get('name') or '').replace('"', '\\"')),
        ]))
    return ';'.join(parts)


def build_compound_doc(cpd):
    doc = dict(cpd)
    for k in DROP_FROM_COMPOUND:
        doc.pop(k, None)
    return doc


def build_reaction_doc(rxn):
    doc = dict(rxn)
    for k in DROP_FROM_REACTION:
        doc.pop(k, None)
    stoich = rxn.get('stoichiometry')
    if isinstance(stoich, list):
        doc['stoichiometry'] = flatten_stoichiometry(stoich)
    return doc


def main():
    compound_docs = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM_ROOT, 'compound_*.json'))):
        with open(path) as fh:
            for cpd in json.load(fh):
                compound_docs.append(build_compound_doc(cpd))
    with open(OUT_COMPOUNDS, 'w') as fh:
        json.dump(compound_docs, fh, indent=None, separators=(',', ':'))

    reaction_docs = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM_ROOT, 'reaction_*.json'))):
        with open(path) as fh:
            for rxn in json.load(fh):
                reaction_docs.append(build_reaction_doc(rxn))
    with open(OUT_REACTIONS, 'w') as fh:
        json.dump(reaction_docs, fh, indent=None, separators=(',', ':'))

    print(f'Wrote {len(compound_docs):,} compounds to {os.path.basename(OUT_COMPOUNDS)}')
    print(f'Wrote {len(reaction_docs):,} reactions to {os.path.basename(OUT_REACTIONS)}')


if __name__ == '__main__':
    main()
