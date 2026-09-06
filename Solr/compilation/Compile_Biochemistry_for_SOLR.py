#!/usr/bin/env python3
"""Compile the biochemistry JSON into two Solr-ready payloads with
nested child documents.

Reads:   ../../Biochemistry/compound_*.json
         ../../Biochemistry/reaction_*.json
Writes:  ./solr_compounds.json
         ./solr_reactions.json

The output is designed for Solr 9 with classic-schema + nested docs
(see ../cores/{compounds,reactions}/schema.xml). Children are emitted
using **named subdocuments** — arrays under keys like `thermodynamics`,
`stoichiometry`, and `pkas` at the parent's top level. This form lets
Solr auto-populate `_nest_path_` (to `/thermodynamics`, etc.) and
`_root_` (parent's unique key) so the `[child]` transformer and
`{!parent}` / `{!child}` block-join query parsers work out of the box.
Each child also carries an explicit `doc_type` so filters can read
readably as `doc_type:thermodynamics` in addition to `_nest_path_:...`.

Denormalized flat fields are computed at compile time so common
browse-page filters (has_structure, has_atom_mapping, atom_count_C,
n_sources_thermodynamics, sources_agree_direction, is_electron_carrier,
n_reactants / n_products) don't require a parent/child join at query
time.
"""
import glob
import json
import os
import re
from collections import Counter

BIOCHEM_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            '..', '..', 'Biochemistry'))
OUT_COMPOUNDS = os.path.join(os.path.dirname(__file__), 'solr_compounds.json')
OUT_REACTIONS = os.path.join(os.path.dirname(__file__), 'solr_reactions.json')

# --- Denormalization helpers -------------------------------------------------

# Electron-carrier anchor compounds (Fdx / Trx / Flx / Grx pairs). Match the
# CONVENTION_A_ANCHORS additions on the `electron-carrier-anchors` branch;
# safe to include here now because the flag is derived from ID membership,
# not from the anchor's numeric ΔfG value.
ELECTRON_CARRIERS = frozenset({
    'cpd27757', 'cpd28082',  # Fdx-ox / Fdx-red   (Fe2R8S2)
    'cpd11621', 'cpd11620',  # Fdx-ox2 / Fdx-red2 (Fe2R4S6)
    'cpd27735', 'cpd28060',  # Trx-ox / Trx-red
    'cpd48855', 'cpd48919',  # Trx-ox-1 / Trx-red-1
    'cpd12207', 'cpd12173',  # Flx-ox / Flx-red
    'cpd15464', 'cpd15465',  # Flx-ox-X / Flx-red-X
    'cpd15480', 'cpd15481',  # Grx-ox / Grx-red
    'cpd27732', 'cpd28057',  # Grx-ox-g / Grx-red-g
})

_ELEMENT_RE = re.compile(r'([A-Z][a-z]?)(\d*)')


def parse_formula_atoms(formula):
    """Return a Counter of element → atom count from a chemical formula.
    Handles two-letter elements and implicit counts of 1.
    Returns an empty Counter for empty / unparseable formulas."""
    counts = Counter()
    if not formula:
        return counts
    for match in _ELEMENT_RE.finditer(formula):
        el = match.group(1)
        if not el:
            continue
        n = int(match.group(2)) if match.group(2) else 1
        counts[el] += n
    return counts


def as_bool(value):
    """Solr's BoolField accepts true/false JSON literals. Convert 0/1
    ints and None to booleans; leave real booleans alone."""
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return bool(value)
    # Accept "1"/"0" / "true"/"false" strings too, defensively.
    if isinstance(value, str):
        s = value.strip().lower()
        if s in ('1', 'true', 't', 'yes'):
            return True
        if s in ('0', 'false', 'f', 'no', ''):
            return False
    return None


def sources_agree_direction(thermo):
    """True when every source's operator is the same non-'?' value.
    Returns None when <2 sources have a non-'?' operator (undefined).

    thermo: dict source → [energy, error, operator]
    """
    ops = []
    for _, val in (thermo or {}).items():
        if not isinstance(val, list) or len(val) < 3:
            continue
        op = val[2]
        if op in ('>', '<', '='):
            ops.append(op)
    if len(ops) < 2:
        return None
    return len(set(ops)) == 1


# --- Compound compilation ----------------------------------------------------

def build_compound_doc(cpd):
    """Return a Solr-ready compound doc (parent + nested children)."""
    doc = {}

    # Flat pass-throughs. Skip fields Solr can't index cleanly (nested
    # dicts we handle separately below; None-valued fields we omit so
    # Solr doesn't have to store a null placeholder).
    _FLAT_STR_KEYS = (
        'id', 'name', 'abbreviation', 'formula', 'source', 'smiles',
        'inchikey', 'linked_compound', 'abstract_compound', 'comprised_of',
        'class',
    )
    _FLAT_NUM_KEYS = ('mass', 'charge', 'deltag', 'deltagerr')
    _FLAT_BOOL_KEYS = ('is_core', 'is_obsolete', 'is_cofactor')
    _FLAT_LIST_KEYS = ('aliases', 'notes', 'pka', 'pkb')

    for k in _FLAT_STR_KEYS:
        v = cpd.get(k)
        if v is not None and v != '':
            doc[k] = v
    for k in _FLAT_NUM_KEYS:
        v = cpd.get(k)
        if v is not None:
            doc[k] = v
    for k in _FLAT_BOOL_KEYS:
        v = as_bool(cpd.get(k))
        if v is not None:
            doc[k] = v
    for k in _FLAT_LIST_KEYS:
        v = cpd.get(k)
        if v:
            if isinstance(v, str):
                v = [v]
            doc[k] = list(v)

    # Explicit parent doc_type — makes block-join queries readable.
    doc['doc_type'] = 'compound'

    # Denormalized derived fields
    atoms = parse_formula_atoms(cpd.get('formula') or '')
    doc['has_structure'] = bool((cpd.get('smiles') or '').strip())
    doc['atom_count_C'] = int(atoms.get('C', 0))
    doc['atom_count_N'] = int(atoms.get('N', 0))
    doc['atom_count_O'] = int(atoms.get('O', 0))
    doc['atom_count_P'] = int(atoms.get('P', 0))
    doc['atom_count_S'] = int(atoms.get('S', 0))
    doc['is_electron_carrier'] = cpd['id'] in ELECTRON_CARRIERS

    thermo = cpd.get('thermodynamics') if isinstance(cpd.get('thermodynamics'), dict) else {}
    doc['n_sources_thermodynamics'] = len(thermo)

    # Named subdocs — Solr auto-populates _nest_path_="/thermodynamics"
    # and _root_=<parent id> for each entry in this array.
    thermo_children = []
    for source, val in sorted(thermo.items()):
        if not isinstance(val, list) or len(val) < 2:
            continue
        child = {
            'id': f"{cpd['id']}::thermo::{source}",
            'doc_type': 'thermodynamics',
            'source_name': source,
        }
        try:
            child['energy'] = float(val[0])
        except (TypeError, ValueError):
            pass
        try:
            child['error'] = float(val[1])
        except (TypeError, ValueError):
            pass
        thermo_children.append(child)
    if thermo_children:
        doc['thermodynamics'] = thermo_children

    pkas = cpd.get('pkas') if isinstance(cpd.get('pkas'), dict) else {}
    pkas_children = []
    for source, val in sorted(pkas.items()):
        if not isinstance(val, dict):
            continue
        child = {
            'id': f"{cpd['id']}::pkas::{source}",
            'doc_type': 'pkas',
            'source_name': source,
        }
        pka = val.get('pKa') or val.get('pka')
        pkb = val.get('pKb') or val.get('pkb')
        if pka is not None:
            child['pka_value'] = [pka] if isinstance(pka, str) else list(pka)
        if pkb is not None:
            child['pkb_value'] = [pkb] if isinstance(pkb, str) else list(pkb)
        pkas_children.append(child)
    if pkas_children:
        doc['pkas'] = pkas_children

    return doc


# --- Reaction compilation ----------------------------------------------------

def build_reaction_doc(rxn):
    """Return a Solr-ready reaction doc (parent + nested children)."""
    doc = {}

    _FLAT_STR_KEYS = (
        'id', 'name', 'abbreviation', 'code', 'equation', 'definition',
        'status', 'source', 'linked_reaction', 'abstract_reaction',
        'reversibility',
    )
    _FLAT_NUM_KEYS = ('deltag', 'deltagerr')
    _FLAT_BOOL_KEYS = ('is_transport', 'is_obsolete')
    _FLAT_LIST_KEYS = ('pathways', 'aliases', 'ec_numbers', 'notes')

    for k in _FLAT_STR_KEYS:
        v = rxn.get(k)
        if v is not None and v != '':
            doc[k] = v
    for k in _FLAT_NUM_KEYS:
        v = rxn.get(k)
        if v is not None:
            doc[k] = v
    for k in _FLAT_BOOL_KEYS:
        v = as_bool(rxn.get(k))
        if v is not None:
            doc[k] = v
    for k in _FLAT_LIST_KEYS:
        v = rxn.get(k)
        if v:
            if isinstance(v, str):
                v = [v]
            doc[k] = list(v)

    doc['doc_type'] = 'reaction'

    # compound_ids on the JSON is a semicolon-delimited string; Solr wants
    # a multi-valued string field, so split it here.
    compound_ids = rxn.get('compound_ids')
    if isinstance(compound_ids, str) and compound_ids:
        doc['compound_ids'] = compound_ids.split(';')
    elif isinstance(compound_ids, list):
        doc['compound_ids'] = compound_ids

    # atom_mapping is a single nested dict on the JSON side:
    #   {"data": [...], "confidence": "clean|salvaged", "has_symmetry_groups": bool}
    # Flatten it out into three sibling fields for Solr indexing so the
    # UI can filter/browse without traversing a nested field name at query
    # time. Denormalized `has_atom_mapping` remains for quick-existence
    # filtering. If the dict is absent (reaction not in the clean set),
    # none of these fields are emitted.
    am = rxn.get('atom_mapping')
    if isinstance(am, dict):
        doc['atom_mapping_data'] = list(am.get('data') or [])
        conf = am.get('confidence')
        if conf:
            doc['atom_mapping_confidence'] = conf
        doc['atom_mapping_has_symmetry_groups'] = bool(am.get('has_symmetry_groups'))
        doc['has_atom_mapping'] = bool(am.get('data'))
    else:
        doc['has_atom_mapping'] = False
    thermo = rxn.get('thermodynamics') if isinstance(rxn.get('thermodynamics'), dict) else {}
    doc['n_sources_thermodynamics'] = len(thermo)
    agree = sources_agree_direction(thermo)
    if agree is not None:
        doc['sources_agree_direction'] = agree

    stoich_entries = rxn.get('stoichiometry') or []
    n_reactants = sum(1 for e in stoich_entries
                      if isinstance(e, dict) and (e.get('coefficient') or 0) < 0)
    n_products = sum(1 for e in stoich_entries
                     if isinstance(e, dict) and (e.get('coefficient') or 0) > 0)
    doc['n_reactants'] = n_reactants
    doc['n_products'] = n_products

    # Named subdocs — thermodynamics per source
    thermo_children = []
    for source, val in sorted(thermo.items()):
        if not isinstance(val, list) or len(val) < 2:
            continue
        child = {
            'id': f"{rxn['id']}::thermo::{source}",
            'doc_type': 'thermodynamics',
            'source_name': source,
        }
        try:
            child['energy'] = float(val[0])
        except (TypeError, ValueError):
            pass
        try:
            child['error'] = float(val[1])
        except (TypeError, ValueError):
            pass
        if len(val) >= 3 and val[2] in ('>', '<', '=', '?'):
            child['operator'] = val[2]
        thermo_children.append(child)
    if thermo_children:
        doc['thermodynamics'] = thermo_children

    # Named subdocs — stoichiometry per participant
    stoich_children = []
    for idx, entry in enumerate(stoich_entries):
        if not isinstance(entry, dict):
            continue
        compound = entry.get('compound')
        if not compound:
            continue
        coef = entry.get('coefficient') or 0
        child = {
            'id': f"{rxn['id']}::stoich::{idx}",
            'doc_type': 'stoichiometry',
            'compound': compound,
            'coefficient': int(coef),
            'is_reactant': coef < 0,
        }
        if entry.get('compartment') is not None:
            child['compartment'] = int(entry['compartment'])
        if entry.get('charge') is not None:
            child['participant_charge'] = int(entry['charge'])
        if entry.get('formula'):
            child['participant_formula'] = entry['formula']
        if entry.get('name'):
            child['participant_name'] = entry['name']
        stoich_children.append(child)
    if stoich_children:
        doc['stoichiometry'] = stoich_children

    return doc


# --- Main --------------------------------------------------------------------

def main():
    # Compounds
    compound_docs = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM_ROOT, 'compound_*.json'))):
        with open(path) as fh:
            for cpd in json.load(fh):
                compound_docs.append(build_compound_doc(cpd))
    with open(OUT_COMPOUNDS, 'w') as fh:
        json.dump(compound_docs, fh, indent=None, separators=(',', ':'))

    # Reactions
    reaction_docs = []
    for path in sorted(glob.glob(os.path.join(BIOCHEM_ROOT, 'reaction_*.json'))):
        with open(path) as fh:
            for rxn in json.load(fh):
                reaction_docs.append(build_reaction_doc(rxn))
    with open(OUT_REACTIONS, 'w') as fh:
        json.dump(reaction_docs, fh, indent=None, separators=(',', ':'))

    def count_children(doc, keys):
        return sum(len(doc.get(k) or []) for k in keys)
    n_c_children = sum(count_children(d, ('thermodynamics', 'pkas')) for d in compound_docs)
    n_r_children = sum(count_children(d, ('thermodynamics', 'stoichiometry')) for d in reaction_docs)
    print(f'Wrote {len(compound_docs):,} compound parents (+{n_c_children:,} child docs) to {OUT_COMPOUNDS}')
    print(f'Wrote {len(reaction_docs):,} reaction parents (+{n_r_children:,} child docs) to {OUT_REACTIONS}')


if __name__ == '__main__':
    main()
