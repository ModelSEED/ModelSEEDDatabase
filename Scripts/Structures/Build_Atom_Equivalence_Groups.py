#!/usr/bin/env python3
"""Build per-compound atom equivalence groups from InChI AuxInfo /E: layer.

Reads:   Biochemistry/Structures/Unique_ModelSEED_Structures.txt (SMILES per cpd)
Writes:  Biochemistry/Structures/AtomMappings/species_equivalence_groups.tsv

Motivation
==========

RDT's atom-mapping output for reactions involving symmetric groups
(CO2's two Os in a decarboxylation, PPi's six equivalent terminal Os,
benzene's six equivalent aromatic Cs) collapses into ambiguous run-on
chains like `A=B=C=D`. Our current row-level filter
(Rebuild_AtomMappings_from_raw.py) splits these into safe pairs but
loses the "these atoms are interchangeable, don't care which specific
pairing" semantics.

InChI's AuxInfo /E: layer already encodes this — atoms that are
equivalent under InChI's canonicalization algorithm. Sebastian Huhn's
MetaCyc-side prototype at
    https://github.com/sebahu/UniversalRDT/tree/main/MetaCyc
uses three shell scripts (create_inchi_equivalent_atoms{,2,3}.sh) to
extract and re-project these groups onto RDT's mapping-row notation.
This script is the Python translation for ModelSEED.

Method
======

For each compound with a SMILES:

  1. Convert SMILES -> InChI + AuxInfo via OpenBabel (`obabel -oinchi -xa`).
     Same tool RDT uses for structure processing, so canonical atom
     numbering agrees with the numbering RDT emits in its mapping rows.
  2. Parse AuxInfo for the /E: layer — a paren+comma nested list like
     `(1,2,3,4,5,6)(8,9)` giving 1-indexed InChI atom numbers grouped
     by equivalence class. Compounds with no /E: layer have no symmetry
     that InChI could detect and are skipped.
  3. Reproject each InChI atom number to RDT's `Element#K` notation:
     enumerate non-H elements alphabetically as they appear in the
     compound formula, and within each element assign sequential K
     starting at 1. This mirrors how InChI's canonical numbering
     orders atoms (C first, then hetero-atoms alphabetically), which
     matches OpenBabel/RDT's `Element#K` scheme.

Output format (tab-separated):
    cpd_id\tequivalence_groups

where equivalence_groups is space-separated set notation:
    cpd00011\t(O#1;O#2)
    cpd00012\t(O#1;O#2;O#3;O#4;O#5;O#6) (P#1;P#2)
    cpd00020\t(O#2;O#3)

Downstream (Rebuild_AtomMappings_from_raw.py) uses this table to
rewrite atom-mapping rows: any reference to an atom in an equivalence
group is replaced with the whole group in set notation.
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
from collections import Counter

# pybel keeps OpenBabel in-process — ~50× faster than one subprocess per
# compound. Requires `openbabel` conda package (bundled with the
# `universalrdt` micromamba env).
try:
    from openbabel import pybel
except ImportError:
    pybel = None

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
IN_STRUCT = os.path.join(REPO, 'Biochemistry', 'Structures',
                         'Unique_ModelSEED_Structures.txt')
OUT_TSV = os.path.join(REPO, 'Biochemistry', 'Structures', 'AtomMappings',
                       'species_equivalence_groups.tsv')

# /E: is nested: outer parens delimit each equivalence class; inner
# comma-list gives 1-indexed InChI atom numbers.  Example: /E:(1,2,3)(5,6)
E_LAYER_RE = re.compile(r'/E:((?:\([\d,]+\))+)')
E_GROUP_RE = re.compile(r'\(([\d,]+)\)')
FORMULA_RE = re.compile(r'([A-Z][a-z]?)(\d*)')


def load_compound_smiles(path):
    """Return {cpd_id: smiles} for each compound with a SMILES row.
    File is 6-col TSV: id, type, aliases, formula, charge, structure."""
    out = {}
    with open(path) as fh:
        header = fh.readline()  # discard
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 6 and parts[1] == 'SMILE':
                cpd = parts[0]
                smiles = parts[5]
                if smiles and cpd not in out:
                    out[cpd] = smiles
    return out


def smiles_to_inchi_auxinfo(smiles):
    """Return (inchi, auxinfo) tuple via pybel (in-process OpenBabel).
    Returns (None, None) on any parse failure."""
    if pybel is None:
        return None, None
    try:
        mol = pybel.readstring('smi', smiles)
        out = mol.write('inchi', opt={'a': True})
    except Exception:
        return None, None
    lines = [l.strip() for l in out.splitlines() if l.strip()]
    inchi = next((l for l in lines if l.startswith('InChI=')), None)
    auxinfo = next((l for l in lines if l.startswith('AuxInfo=')), None)
    return inchi, auxinfo


def parse_equivalence_groups(auxinfo):
    """Return [[int, ...], ...] of 1-indexed InChI atom numbers, one
    list per equivalence class. Empty list if no /E: layer."""
    if not auxinfo:
        return []
    m = E_LAYER_RE.search(auxinfo)
    if not m:
        return []
    groups = []
    for g in E_GROUP_RE.finditer(m.group(1)):
        groups.append([int(x) for x in g.group(1).split(',') if x])
    return groups


def inchi_atom_to_element_k(formula):
    """Return [Element#K string, ...] indexed by 1-based InChI atom number.
    Order follows the alphabetical-formula-enumeration heuristic Sebastian
    uses: parse formula, sort element symbols alphabetically (skipping H),
    then within each element assign K=1..N. Returns None-padded so index 0
    is unused (InChI atom numbers are 1-based).

    Empty list if formula is empty or unparseable."""
    if not formula:
        return []
    element_counts = Counter()
    for match in FORMULA_RE.finditer(formula):
        el = match.group(1)
        if el:
            n = int(match.group(2)) if match.group(2) else 1
            element_counts[el] += n
    result = [None]  # index 0 unused
    for el in sorted(k for k in element_counts if k != 'H'):
        for k in range(1, element_counts[el] + 1):
            result.append('{}#{}'.format(el, k))
    return result


def process_compound(cpd_id, smiles, formula_from_inchi):
    """Return equivalence-groups string (space-separated set notation)
    for one compound, or None if the compound has no /E: layer or the
    conversion failed."""
    inchi, auxinfo = smiles_to_inchi_auxinfo(smiles)
    if not auxinfo:
        return None
    groups = parse_equivalence_groups(auxinfo)
    if not groups:
        return None
    # InChI formula is the /... segment right after "InChI=1S/"
    if formula_from_inchi is None:
        m = re.match(r'InChI=1S?/([^/]+)/', inchi or '')
        formula_from_inchi = m.group(1) if m else ''
    atom_labels = inchi_atom_to_element_k(formula_from_inchi)
    if not atom_labels:
        return None
    set_strings = []
    for grp in groups:
        try:
            labels = [atom_labels[i] for i in grp]
        except IndexError:
            continue  # skip malformed groups pointing past our atom count
        if all(labels):
            set_strings.append('(' + ';'.join(labels) + ')')
    if not set_strings:
        return None
    return ' '.join(set_strings)


def main():
    if not os.path.isfile(IN_STRUCT):
        sys.exit('Missing input: ' + IN_STRUCT)

    smiles_by_cpd = load_compound_smiles(IN_STRUCT)
    n_total = len(smiles_by_cpd)
    print('Loaded {:,} compounds with SMILES from {}'.format(
        n_total, os.path.basename(IN_STRUCT)), file=sys.stderr)

    n_processed = 0
    n_with_groups = 0
    n_conversion_failed = 0
    results = {}

    for i, (cpd, smi) in enumerate(sorted(smiles_by_cpd.items())):
        n_processed += 1
        inchi, aux = smiles_to_inchi_auxinfo(smi)
        if inchi is None:
            n_conversion_failed += 1
            continue
        m = re.match(r'InChI=1S?/([^/]+)/', inchi)
        formula = m.group(1) if m else ''
        groups_str = process_compound(cpd, smi, formula)
        if groups_str:
            results[cpd] = groups_str
            n_with_groups += 1
        if (i + 1) % 500 == 0:
            print('  processed {:,}/{:,} ({:,} with groups)'.format(
                i + 1, n_total, n_with_groups), file=sys.stderr)

    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, 'w') as fh:
        fh.write('cpd_id\tequivalence_groups\n')
        for cpd in sorted(results):
            fh.write('{}\t{}\n'.format(cpd, results[cpd]))

    print('', file=sys.stderr)
    print('Compounds processed:            {:,}'.format(n_processed), file=sys.stderr)
    print('Compounds obabel failed on:     {:,}'.format(n_conversion_failed), file=sys.stderr)
    print('Compounds with /E: groups:      {:,} ({:.1f}%)'.format(
        n_with_groups, 100 * n_with_groups / max(1, n_processed)), file=sys.stderr)
    print('Wrote {} equivalence-group rows to {}'.format(
        len(results), os.path.basename(OUT_TSV)), file=sys.stderr)


if __name__ == '__main__':
    main()
