#!/usr/bin/env python3
"""Compile the structures core payload: one doc per compound_id with
SMILES, InChI, InChIKey, and a pre-rendered 2D-structure SVG.

Reads:   ../../Biochemistry/compound_*.json          (compound id list)
         ../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt
                                                     (SMILES / InChI / InChIKey per cpd)
Writes:  ./solr_structures.json

Purpose
=======

The compound / reaction cores focus on metabolic-model semantics —
formulas, thermodynamics, atom mappings. The `structures` core is a
dedicated place for structure-representation data, keyed by
compound_id, so a compound-detail page's structure viewer does one
GET on /solr/structures/select?q=id:cpdXXXXX and gets back everything
it needs without wading through the compounds core's per-compound
metadata.

Design
======

One document per compound_id (all ~45.7K compounds, whether they have
a structure or not). Fields are optional — compounds without SMILES
just have the `smiles` field absent. Existence filters use
`smiles:[* TO *]` and negation `-smiles:[* TO *]`.

The `svg` field carries a pre-rendered RDKit MolDraw2DSVG at 300x300.
Rendering moves the CPU cost from UI-serve-time to compile-time
(once, at container-build). Skipped when:
  - the compound has no SMILES
  - the SMILES contains `*` (wildcard R-groups render as literal
    `*` atoms — ugly and misleading)
  - the rendered SVG exceeds 100 KB (rare — <0.5% of compounds,
    mostly polymers whose depiction wouldn't be useful at
    thumbnail size anyway)
Consumers can regenerate SVG on the fly for the skipped compounds.

Same JSON is posted to both `structures_staging` and (eventually)
bare `structures` cores — the structures schema doesn't have a
legacy variant because the production UI doesn't have this core
today; nothing to be backward-compatible with.
"""
import glob
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

# RDKit is chatty about non-fatal parse warnings; quiet it down for
# the batch run.
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.info')

BIOCHEM_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                            '..', '..', 'Biochemistry'))
UNIQUE_STRUCT = os.path.join(BIOCHEM_ROOT, 'Structures',
                             'Unique_ModelSEED_Structures.txt')
OUT = os.path.join(os.path.dirname(__file__), 'solr_structures.json')

SVG_WIDTH = 300
SVG_HEIGHT = 300
SVG_MAX_BYTES = 100 * 1024   # 100 KB cap; skip larger renderings


def load_structures(path):
    """Return {cpd_id: {'smiles': str, 'inchi': str, 'inchikey': str}}.
    Each compound has whatever representations are present in the
    Unique_ModelSEED_Structures.txt file — some have all three, some
    only SMILES + InChIKey, some none."""
    out = {}
    with open(path) as fh:
        fh.readline()  # discard header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            cpd, stype, _aliases, _formula, _charge, structure = parts[:6]
            key = {'SMILE': 'smiles', 'InChI': 'inchi',
                   'InChIKey': 'inchikey'}.get(stype)
            if key and structure:
                out.setdefault(cpd, {}).setdefault(key, structure)
    return out


def load_compound_ids(biochem_root):
    """Return a sorted list of ALL compound_ids in the DB, whether or
    not they have a stored structure. Ensures the structures core has
    a doc for every compound so a UI query never returns 'not found'
    for a known compound."""
    ids = set()
    for path in sorted(glob.glob(os.path.join(biochem_root, 'compound_*.json'))):
        with open(path) as fh:
            for c in json.load(fh):
                if c.get('id'):
                    ids.add(c['id'])
    return sorted(ids)


def render_svg(smiles):
    """Return the SVG text for a SMILES string, or None if the
    compound is unrenderable (parse failure, wildcard atoms, or
    rendering exceeds SVG_MAX_BYTES)."""
    if not smiles or '*' in smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        return None
    if mol is None:
        return None
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(SVG_WIDTH, SVG_HEIGHT)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
    except Exception:
        return None
    if len(svg) > SVG_MAX_BYTES:
        return None
    return svg


def main():
    print(f'Loading structures table: {os.path.basename(UNIQUE_STRUCT)}', file=sys.stderr)
    structures = load_structures(UNIQUE_STRUCT)
    print(f'  {len(structures):,} compounds with at least one structure representation',
          file=sys.stderr)

    print('Loading compound id list from compound_*.json ...', file=sys.stderr)
    all_ids = load_compound_ids(BIOCHEM_ROOT)
    print(f'  {len(all_ids):,} total compound ids', file=sys.stderr)

    docs = []
    n_svg = 0
    n_svg_skipped_wildcard = 0
    n_svg_skipped_oversize = 0
    n_svg_skipped_no_smiles = 0
    n_svg_skipped_parse = 0
    for i, cpd_id in enumerate(all_ids):
        doc = {'id': cpd_id}
        s = structures.get(cpd_id, {})
        if s.get('smiles'):
            doc['smiles'] = s['smiles']
        if s.get('inchi'):
            doc['inchi'] = s['inchi']
        if s.get('inchikey'):
            doc['inchikey'] = s['inchikey']

        smi = s.get('smiles')
        if not smi:
            n_svg_skipped_no_smiles += 1
        elif '*' in smi:
            n_svg_skipped_wildcard += 1
        else:
            svg = render_svg(smi)
            if svg is None:
                # Could be parse failure or oversize — classify
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    n_svg_skipped_parse += 1
                else:
                    n_svg_skipped_oversize += 1
            else:
                doc['svg'] = svg
                n_svg += 1

        docs.append(doc)
        if (i + 1) % 5000 == 0:
            print(f'  processed {i+1:,}/{len(all_ids):,} ({n_svg:,} SVGs rendered)',
                  file=sys.stderr)

    with open(OUT, 'w') as fh:
        json.dump(docs, fh, indent=None, separators=(',', ':'))

    total_svg_bytes = sum(len(d.get('svg') or '') for d in docs)
    print(f'\nWrote {len(docs):,} docs to {os.path.basename(OUT)}', file=sys.stderr)
    print(f'  with SVG:           {n_svg:>7,}', file=sys.stderr)
    print(f'  no SMILES:          {n_svg_skipped_no_smiles:>7,}', file=sys.stderr)
    print(f'  wildcard SMILES:    {n_svg_skipped_wildcard:>7,}', file=sys.stderr)
    print(f'  SVG oversize (>{SVG_MAX_BYTES//1024} KB): {n_svg_skipped_oversize:>7,}', file=sys.stderr)
    print(f'  SMILES parse fail:  {n_svg_skipped_parse:>7,}', file=sys.stderr)
    print(f'  total SVG bytes:    {total_svg_bytes/(1024*1024):.1f} MB', file=sys.stderr)


if __name__ == '__main__':
    main()
