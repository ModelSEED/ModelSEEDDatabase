#!/usr/bin/env python
"""
Refresh formula and charge columns in the per-source structure files
by re-parsing each structure with RDKit (preferred) or OpenBabel.

Reads/writes (in-place) the post-A1 layout:

  Biochemistry/Structures/<source>/inchi.tsv
                                  /smiles.tsv
                                  /protonations/<tool>_<ver>_ph<n>.tsv

inchikey.tsv is left alone (an InChIKey is a hash, no formula/charge to
derive). The protonations file holds rows for InChI, SMILE, and
InChIKey — only the first two have their formula/charge refreshed; the
InChIKey rows are passed through unchanged.

You should re-run this script whenever a structure changes (new KEGG
release, curator edit, new Marvin protonation), or when RDKit/OpenBabel
is upgraded — different parser versions can produce slightly different
formulas/charges. At time of writing we use RDKit 2022.03.5 and
OpenBabel 3.1.1 (same as in the previous script generation).

Outputs diagnostic reports to Biochemistry/Structures/_reports/.
"""
import csv
import glob
import os
import re
import sys

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds  # noqa: E402

from openbabel import pybel  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402
from rdkit import RDLogger  # noqa: E402

RDLogger.logger().setLevel(RDLogger.ERROR)

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
STRUCT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..',
                                            'Biochemistry', 'Structures'))
REPORT_DIR  = os.path.join(STRUCT_ROOT, '_reports')

SOURCES = ['KEGG', 'MetaCyc', 'ChEBI', 'Rhea']


def parse_structure(struct_type, structure):
    """Return (formula, charge, mol_source) or (None, None, None) if neither
    parser succeeds. mol_source is 'RDKit' or 'OpenBabel' depending on
    which one produced the result (RDKit preferred).
    """
    mol_rdkit  = None
    mol_obabel = None
    try:
        if struct_type == 'InChI':
            mol_rdkit  = AllChem.MolFromInchi(structure)
            mol_obabel = pybel.readstring('inchi', structure)
        elif struct_type == 'SMILE':
            mol_rdkit  = AllChem.MolFromSmiles(structure)
            mol_obabel = pybel.readstring('smiles', structure)
    except Exception:
        pass

    if mol_rdkit is None and mol_obabel is None:
        return None, None, None

    if mol_rdkit is not None:
        formula = AllChem.CalcMolFormula(mol_rdkit)
        charge  = AllChem.GetFormalCharge(mol_rdkit)
        mol_src = 'RDKit'
        m = re.search(r'([-+]\d?)$', formula)
        if m:
            formula = formula.replace(m.group(), '')
    else:
        formula = mol_obabel.formula
        charge  = mol_obabel.charge
        mol_src = 'OpenBabel'
        m = re.search(r'([-+]+)$', formula)
        if m:
            formula = formula.replace(m.group(), '')

    formula = re.sub(r'\*', 'R', formula)
    formula = Compounds.mergeFormula(formula)[0]
    return formula, str(charge), mol_src


def refresh_file(path, struct_type_for_all=None, type_column='type',
                 structure_column='structure', report_resolved=None,
                 report_unresolved=None, source_label=''):
    """Refresh formula/charge columns in path, in-place.

    - If struct_type_for_all is set (e.g. 'InChI' or 'SMILE'), every row's
      structure_column is parsed as that type. Used for inchi.tsv and
      smiles.tsv where each file is single-type.
    - Otherwise, the row's type_column tells us how to parse. Used for
      protonations/*.tsv where each row carries its own type.
    """
    if not os.path.isfile(path):
        return

    with open(path) as fh:
        reader     = csv.DictReader(fh, dialect='excel-tab')
        fieldnames = reader.fieldnames
        rows       = list(reader)

    for row in rows:
        struct = row.get(structure_column, '')
        stype  = struct_type_for_all if struct_type_for_all else row.get(type_column)
        ext_id = row.get('external_id') or row.get('ID') or ''
        if stype == 'InChIKey' or not struct or not stype:
            continue

        formula, charge, mol_src = parse_structure(stype, struct)
        if formula is None:
            if report_unresolved is not None:
                report_unresolved.write('\t'.join([source_label, ext_id, stype, struct]) + '\n')
            continue

        row['formula'] = formula
        row['charge']  = charge
        if report_resolved is not None:
            report_resolved.write('\t'.join([source_label, ext_id, stype, struct, formula, charge, mol_src]) + '\n')

    # Write back in the same plain-tab/LF format Migrate_To_New_Layout
    # used. csv.DictWriter('excel-tab') writes \r\n which would diverge
    # from the rest of the repo, so format manually.
    with open(path, 'w') as fh:
        fh.write('\t'.join(fieldnames) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row.get(f, '') or '') for f in fieldnames) + '\n')


def main():
    os.makedirs(REPORT_DIR, exist_ok=True)
    resolved   = open(os.path.join(REPORT_DIR, 'Resolved_Structures.txt'),   'w')
    unresolved = open(os.path.join(REPORT_DIR, 'Unresolved_Structures.txt'), 'w')
    try:
        for source in SOURCES:
            src_dir = os.path.join(STRUCT_ROOT, source)
            print(f'Refreshing {source}...')
            refresh_file(os.path.join(src_dir, 'inchi.tsv'),  struct_type_for_all='InChI',
                         report_resolved=resolved, report_unresolved=unresolved,
                         source_label=source)
            refresh_file(os.path.join(src_dir, 'smiles.tsv'), struct_type_for_all='SMILE',
                         report_resolved=resolved, report_unresolved=unresolved,
                         source_label=source)

            proto_dir = os.path.join(src_dir, 'protonations')
            if os.path.isdir(proto_dir):
                for proto_file in sorted(glob.glob(os.path.join(proto_dir, '*.tsv'))):
                    refresh_file(proto_file,
                                 report_resolved=resolved, report_unresolved=unresolved,
                                 source_label=source)
    finally:
        resolved.close()
        unresolved.close()


if __name__ == '__main__':
    main()
