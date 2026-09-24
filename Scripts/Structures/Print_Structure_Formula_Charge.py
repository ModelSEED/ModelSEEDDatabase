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

INCHI ROWS ARE READ FROM THEIR OWN LAYERS, not from a parsed molecule. A
standard InChI declares its formula (the formula layer, with component
multipliers), its protonation (/p) and its charges (/q) outright; the
molecule a parser builds from it can only agree with those or be wrong.
Two ways of being wrong were found in this repository's files:

  * RDKit rejects hypervalent halogen oxides (chlorate, chlorite, iodate,
    `InChI=1S/ClO3/c2-1(3)4/q-1`: "explicit valence for Cl, 6"), the
    OpenBabel fallback warns "Charge(s): Do not match" and returns
    charge 0. Nine rows stored a neutral radical where the string
    declares an anion.
  * On multi-component strings with both /q and /p layers -- the Mg
    porphyrins, `.../q-1;+2/p-1` -- RDKit keeps the pre-/p hydrogen
    count and its charge, so 49 rows stored one H and one charge unit
    too many. The string declares a neutral chlorophyll; RDKit shipped
    a cation.

Both are parser defects, not chemistry, and both propagate into every
protonation bundle refreshed through this function. inchi_layers() below
is the fix: for a standard InChI it returns what the identifier itself
says, deterministically and independently of any parser version, which
also retires the version-drift caveat above for InChI rows. Over the
44,297 standard InChIs in the four source files it agrees with the
parsed result on 44,239 and differs on exactly the 58 above. The
mol-parser path remains for SMILES, and as the fallback for a
non-standard InChI or one with a /f or /r layer.

Outputs diagnostic reports to Biochemistry/Structures/_reports/.
"""

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _parser = _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter)
    _parser.add_argument(
        "--types", nargs="*", choices=("InChI", "SMILE"), default=None,
        help="refresh only rows of these structure types (default: both)")
    _ARGS = _parser.parse_args()


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


_INCHI_ELEMENT = re.compile(r'([A-Z][a-z]?)(\d*)')


def inchi_layers(inchi):
    """(formula, charge) that a standard InChI declares in its own layers.

    formula: the formula layer summed over components ('3C6H5.ClH.Sn' is
             three phenyls, HCl and tin), with the /p proton layer applied
             to the hydrogen count -- /p-1 is one proton removed from the
             structure the formula layer describes.
    charge:  the sum of the /q layer (one entry per component instance,
             'n*' multipliers honoured) plus the /p layer.

    Returns (None, None) for anything that is not a standard InChI, for a
    string carrying a fixed-H (/f) or reconnected (/r) layer -- those are
    non-standard and change the formula -- and for a /p that would take
    the hydrogen count below zero. The caller falls back to a parser.
    """
    if not inchi or not inchi.startswith('InChI=1S/'):
        return None, None
    layers = inchi[len('InChI=1S/'):].split('/')
    if any(l.startswith(('f', 'r')) for l in layers[1:]):
        return None, None
    counts = {}
    for component in layers[0].split('.'):
        m = re.match(r'^(\d*)(.*)$', component)
        mult = int(m.group(1)) if m.group(1) else 1
        for el, n in _INCHI_ELEMENT.findall(m.group(2)):
            counts[el] = counts.get(el, 0) + mult * (int(n) if n else 1)
    protons = charge = 0
    for layer in layers[1:]:
        if layer.startswith('p'):
            protons += int(layer[1:])
        elif layer.startswith('q'):
            for tok in layer[1:].split(';'):
                tok = tok.strip()
                if tok:
                    mult, _, val = tok.rpartition('*')
                    charge += (int(mult) if mult else 1) * int(val)
    counts['H'] = counts.get('H', 0) + protons
    if counts['H'] < 0:
        return None, None
    formula = ''.join(f'{el}{n if n != 1 else ""}' for el, n in counts.items() if n > 0)
    return Compounds.mergeFormula(formula)[0], charge + protons


def parse_structure(struct_type, structure):
    """Return (formula, charge, mol_source) or (None, None, None) if nothing
    can read the structure. mol_source is 'InChI-layers' for a standard
    InChI (see inchi_layers), else 'RDKit' or 'OpenBabel' depending on
    which parser produced the result (RDKit preferred).
    """
    if struct_type == 'InChI':
        formula, charge = inchi_layers(structure)
        if formula is not None:
            return formula, str(charge), 'InChI-layers'

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
    # INVARIANT: a structure carrying wildcard atoms gets an R in its formula.
    # RDKit renders each dummy atom as `*` (rewritten to R just above);
    # OpenBabel omits them entirely, so a SMILES that fell through to it came
    # back without the R the run script had put there -- 8 rows of the 26.1
    # MetaCyc bundle lost their R on a refresh. Count the wildcards in the
    # string itself, which neither parser can lose.
    if struct_type == 'SMILE' and not re.search(r'R(?![a-z])', formula):
        wildcards = structure.count('*')
        if wildcards:
            formula += 'R' if wildcards == 1 else f'R{wildcards}'
    formula = Compounds.mergeFormula(formula)[0]
    return formula, str(charge), mol_src


def refresh_file(path, struct_type_for_all=None, type_column='type',
                 structure_column='structure', report_resolved=None,
                 report_unresolved=None, source_label='', only_types=None):
    """Refresh formula/charge columns in path, in-place.

    - If struct_type_for_all is set (e.g. 'InChI' or 'SMILE'), every row's
      structure_column is parsed as that type. Used for inchi.tsv and
      smiles.tsv where each file is single-type.
    - Otherwise, the row's type_column tells us how to parse. Used for
      protonations/*.tsv where each row carries its own type.
    - only_types restricts the refresh to those structure types; rows of
      any other type are written back untouched.
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
        if only_types and stype not in only_types:
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
    only = getattr(sys.modules[__name__], '_ARGS', None)
    only = only.types if only is not None else None
    os.makedirs(REPORT_DIR, exist_ok=True)
    resolved   = open(os.path.join(REPORT_DIR, 'Resolved_Structures.txt'),   'w')
    unresolved = open(os.path.join(REPORT_DIR, 'Unresolved_Structures.txt'), 'w')
    try:
        for source in SOURCES:
            src_dir = os.path.join(STRUCT_ROOT, source)
            print(f'Refreshing {source}...')
            # The source files name their structure column after the
            # representation -- `inchi`, `smiles` -- not `structure`. With the
            # default column name every source row was skipped and only the
            # bundles were ever refreshed, so the "re-run whenever a structure
            # changes" advice above has not applied to these files since the
            # layout migration. Found because a refresh that should have moved
            # 58 inchi.tsv rows moved none.
            if not only or 'InChI' in only:
                refresh_file(os.path.join(src_dir, 'inchi.tsv'),  struct_type_for_all='InChI',
                             structure_column='inchi',
                             report_resolved=resolved, report_unresolved=unresolved,
                             source_label=source)
            if not only or 'SMILE' in only:
                refresh_file(os.path.join(src_dir, 'smiles.tsv'), struct_type_for_all='SMILE',
                             structure_column='smiles',
                             report_resolved=resolved, report_unresolved=unresolved,
                             source_label=source)

            proto_dir = os.path.join(src_dir, 'protonations')
            if os.path.isdir(proto_dir):
                for proto_file in sorted(glob.glob(os.path.join(proto_dir, '*.tsv'))):
                    refresh_file(proto_file,
                                 report_resolved=resolved, report_unresolved=unresolved,
                                 source_label=source, only_types=only)
    finally:
        resolved.close()
        unresolved.close()


if __name__ == '__main__':
    main()
