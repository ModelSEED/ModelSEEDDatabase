#!/usr/bin/env python
"""
Migrate Biochemistry/Structures/<DB>/ from the legacy layout to the
target layout declared in sources.yaml:

  <DB>/
    inchi.tsv               external_id, inchi, formula, charge
    smiles.tsv              external_id, smiles, formula, charge
    inchikey.tsv            external_id, inchikey
    protonations/
      <tool>_<ver>_ph<n>.tsv  external_id, type, structure, formula,
                              charge, tool, tool_version, ph, generated_on
    pkas/
      <tool>_<ver>.tsv      external_id, kind, value, tool, tool_version

This is the C1 migration step: GENERATE the new layout from the
existing files. Old *_Strings.txt and *_Formulas_Charges.txt files are
left in place so the existing pipeline continues to function. C2
switches readers; C3 retires the old files.

Idempotent: re-running produces byte-identical output. Lossless: all
data in the legacy files is preserved (re-derivable both ways).
"""
import os
import sys

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
STRUCT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..',
                                            'Biochemistry', 'Structures'))

# Snapshot dates recovered from `git log -1 --format=%cI` on the
# Charged structure files. Tool version inferred from
# Curation/ignores/ignored_structures_marvin234_winter2023.txt.
SOURCES = {
    'KEGG':    {'snapshot': '2024-01-23', 'tool': 'Marvin', 'tool_version': '23.4', 'ph': 7},
    'MetaCyc': {'snapshot': '2024-01-23', 'tool': 'Marvin', 'tool_version': '23.4', 'ph': 7},
    'ChEBI':   {'snapshot': '2024-01-23', 'tool': 'Marvin', 'tool_version': '23.4', 'ph': 7},
    'Rhea':    {'snapshot': '2023-03-31', 'tool': 'Marvin', 'tool_version': '23.4', 'ph': 7},
}


def read_strings_file(path):
    """*_Strings.txt: tab-separated (id, structure)."""
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0]:
                out[parts[0]] = parts[1]
    return out


def read_formula_charges_file(path):
    """*_Formulas_Charges.txt: tab-separated (id, formula, charge)."""
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[0]:
                out[parts[0]] = (parts[1], parts[2])
    return out


def read_pka_file(path):
    """pKa_Strings.txt: tab-separated (id, kind, value)."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                rows.append((parts[0], parts[1], parts[2]))
    return rows


def write_tsv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for row in rows:
            fh.write('\t'.join('' if v is None else str(v) for v in row) + '\n')


def migrate_source(source, meta):
    src_dir = os.path.join(STRUCT_ROOT, source)
    print(f'  {source}:')

    inchi_orig = read_strings_file(os.path.join(src_dir, 'InChI_OriginalStrings.txt'))
    smile_orig = read_strings_file(os.path.join(src_dir, 'SMILE_OriginalStrings.txt'))
    ikey_orig  = read_strings_file(os.path.join(src_dir, 'InChIKey_OriginalStrings.txt'))
    inchi_chrg = read_strings_file(os.path.join(src_dir, 'InChI_ChargedStrings.txt'))
    smile_chrg = read_strings_file(os.path.join(src_dir, 'SMILE_ChargedStrings.txt'))
    ikey_chrg  = read_strings_file(os.path.join(src_dir, 'InChIKey_ChargedStrings.txt'))

    inchi_orig_fc = read_formula_charges_file(os.path.join(src_dir, 'InChI_Original_Formulas_Charges.txt'))
    smile_orig_fc = read_formula_charges_file(os.path.join(src_dir, 'SMILE_Original_Formulas_Charges.txt'))
    inchi_chrg_fc = read_formula_charges_file(os.path.join(src_dir, 'InChI_Charged_Formulas_Charges.txt'))
    smile_chrg_fc = read_formula_charges_file(os.path.join(src_dir, 'SMILE_Charged_Formulas_Charges.txt'))

    pka_rows = read_pka_file(os.path.join(src_dir, 'pKa_Strings.txt'))

    # ---- inchi.tsv / smiles.tsv / inchikey.tsv (source-as-downloaded) ----
    rows = []
    for ext_id in sorted(inchi_orig):
        f, c = inchi_orig_fc.get(ext_id, ('', ''))
        rows.append((ext_id, inchi_orig[ext_id], f, c))
    write_tsv(os.path.join(src_dir, 'inchi.tsv'),
              ['external_id', 'inchi', 'formula', 'charge'], rows)
    print(f'    inchi.tsv:    {len(rows)} rows')

    rows = []
    for ext_id in sorted(smile_orig):
        f, c = smile_orig_fc.get(ext_id, ('', ''))
        rows.append((ext_id, smile_orig[ext_id], f, c))
    write_tsv(os.path.join(src_dir, 'smiles.tsv'),
              ['external_id', 'smiles', 'formula', 'charge'], rows)
    print(f'    smiles.tsv:   {len(rows)} rows')

    rows = [(ext_id, ikey_orig[ext_id]) for ext_id in sorted(ikey_orig)]
    write_tsv(os.path.join(src_dir, 'inchikey.tsv'),
              ['external_id', 'inchikey'], rows)
    print(f'    inchikey.tsv: {len(rows)} rows')

    # ---- protonations/<tool>_<ver>_ph<n>.tsv ----
    proto_name = f"{meta['tool'].lower()}_{meta['tool_version']}_ph{meta['ph']}.tsv"
    proto_path = os.path.join(src_dir, 'protonations', proto_name)

    rows = []
    all_charged = sorted(set(inchi_chrg) | set(smile_chrg) | set(ikey_chrg))
    for ext_id in all_charged:
        if ext_id in inchi_chrg:
            f, c = inchi_chrg_fc.get(ext_id, ('', ''))
            rows.append((ext_id, 'InChI', inchi_chrg[ext_id], f, c,
                         meta['tool'], meta['tool_version'], meta['ph'], meta['snapshot']))
        if ext_id in smile_chrg:
            f, c = smile_chrg_fc.get(ext_id, ('', ''))
            rows.append((ext_id, 'SMILE', smile_chrg[ext_id], f, c,
                         meta['tool'], meta['tool_version'], meta['ph'], meta['snapshot']))
        if ext_id in ikey_chrg:
            rows.append((ext_id, 'InChIKey', ikey_chrg[ext_id], '', '',
                         meta['tool'], meta['tool_version'], meta['ph'], meta['snapshot']))
    write_tsv(proto_path,
              ['external_id', 'type', 'structure', 'formula', 'charge',
               'tool', 'tool_version', 'ph', 'generated_on'], rows)
    print(f'    protonations/{proto_name}: {len(rows)} rows')

    # ---- pkas/<tool>_<ver>.tsv ----
    pka_name = f"{meta['tool'].lower()}_{meta['tool_version']}.tsv"
    pka_path = os.path.join(src_dir, 'pkas', pka_name)
    pka_out  = [(ext, kind, val, meta['tool'], meta['tool_version']) for ext, kind, val in pka_rows]
    write_tsv(pka_path,
              ['external_id', 'kind', 'value', 'tool', 'tool_version'], pka_out)
    print(f'    pkas/{pka_name}: {len(pka_out)} rows')


def main():
    print('Generating new layout files...')
    for source, meta in SOURCES.items():
        migrate_source(source, meta)
    print('Done.')


if __name__ == '__main__':
    main()
