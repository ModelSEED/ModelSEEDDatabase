#!/usr/bin/env python
"""
Preview_Structure_Update.py <cpd_id> <new_inchi>

Simulates what would happen if a compound's primary structure were updated
to <new_inchi>, and prints the diff against the current TSV state without
modifying any files.

Reports:
  formula, charge, inchikey, smiles  — derived from RDKit on the new InChI
                                       (matches the production normalization
                                       in Print_Structure_Formula_Charge.py)
  ΔG (eQuilibrator)                  — looked up from the new InChIKey's
                                       first two layers via MetaNetX
  ΔG (Group Contribution)            — NOT computed (would require MFAToolkit
                                       re-run on a new mol file); the script
                                       reports the existing GC value
  pka, pkb                           — informational; pKa attribution is
                                       alias-based and does not change with
                                       structure updates

Requires RDKit (already in the msd-env conda environment used by
Print_Structure_Formula_Charge.py).
"""
import os
import re
import sys

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
BIOCHEM_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))

sys.path.append(os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python'))


def die(msg, code=1):
    print(msg, file=sys.stderr)
    sys.exit(code)


def usage():
    die(__doc__.strip(), code=1)


def diff_tag(a, b):
    return ' [CHANGE]' if str(a) != str(b) else ''


def main():
    if len(sys.argv) != 3:
        usage()
    cpd_id, new_inchi = sys.argv[1], sys.argv[2]

    try:
        from rdkit.Chem import AllChem, inchi
        from rdkit import RDLogger
        RDLogger.logger().setLevel(RDLogger.ERROR)
    except ImportError:
        die("RDKit is required. Activate the msd-env conda environment "
            "(see Scripts/Structures/README.md).")

    from BiochemPy import Compounds

    helper    = Compounds()
    compounds = helper.loadCompounds()
    if cpd_id not in compounds:
        die(f"ERROR: {cpd_id} not in compound dataset", code=2)

    cpd      = compounds[cpd_id]
    cur_thermo = cpd.get('thermodynamics') if isinstance(cpd.get('thermodynamics'), dict) else {}

    print(f"=== {cpd_id}  {cpd.get('name','')!r} ===\n")
    print(f"Current TSV state:")
    print(f"  formula  : {cpd.get('formula')}")
    print(f"  charge   : {cpd.get('charge')}")
    print(f"  inchikey : {cpd.get('inchikey')}")
    print(f"  smiles   : {cpd.get('smiles')}")
    print(f"  pka      : {cpd.get('pka')}")
    print(f"  pkb      : {cpd.get('pkb')}")
    print(f"  deltag   : {cpd.get('deltag')}    deltagerr: {cpd.get('deltagerr')}")
    print(f"  thermodynamics: {cur_thermo}")
    print(f"  notes    : {cpd.get('notes')}")

    # ---- Parse new InChI ------------------------------------------------
    mol = AllChem.MolFromInchi(new_inchi)
    if mol is None:
        die(f"ERROR: RDKit could not parse InChI: {new_inchi}", code=3)

    raw_formula = AllChem.CalcMolFormula(mol)
    new_charge  = AllChem.GetFormalCharge(mol)
    # Strip trailing charge sign from formula (matches Print_Structure_Formula_Charge.py)
    m = re.search(r'([-+]\d?)$', raw_formula)
    if m:
        raw_formula = raw_formula.replace(m.group(), '')
    # Normalize R-groups, Hill sort, merge fragments
    new_formula = re.sub(r'\*', 'R', raw_formula)
    new_formula = Compounds.mergeFormula(new_formula)[0]

    new_inchikey = inchi.InchiToInchiKey(new_inchi)
    new_smiles   = AllChem.MolToSmiles(mol)

    print(f"\nProposed (from new InChI):")
    print(f"  formula  : {new_formula}{diff_tag(new_formula, cpd.get('formula'))}")
    print(f"  charge   : {new_charge}{diff_tag(new_charge, cpd.get('charge'))}")
    print(f"  inchikey : {new_inchikey}{diff_tag(new_inchikey, cpd.get('inchikey'))}")
    print(f"  smiles   : {new_smiles}  (informational — production keeps a separate SMILE pick)")

    # ---- eQuilibrator ΔG lookup ----------------------------------------
    mnx_path = os.path.join(BIOCHEM_ROOT, 'Structures', 'MetaNetX',
                            'Structures_in_ModelSEED_and_eQuilibrator.txt')
    eq_path  = os.path.join(BIOCHEM_ROOT, 'Thermodynamics', 'eQuilibrator',
                            'MetaNetX_Compound_Energies.tbl')

    mnx_by_ikey = {}
    with open(mnx_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                depro = '-'.join(parts[1].split('-')[:2])
                mnx_by_ikey.setdefault(depro, []).append(parts[0])

    eq_energies = {}
    with open(eq_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            if 'energy' in parts[1] or parts[1] == 'nan':
                continue
            try:
                eq_energies[parts[0]] = float(parts[1])
            except ValueError:
                continue

    def eq_lookup(ikey):
        if not ikey:
            return None
        depro = '-'.join(ikey.split('-')[:2])
        if depro not in mnx_by_ikey:
            return None
        best = None
        for mnx in mnx_by_ikey[depro]:
            if mnx in eq_energies:
                if best is None or eq_energies[mnx] < best[1]:
                    best = (mnx, eq_energies[mnx])
        return best

    cur_eq = eq_lookup(cpd.get('inchikey'))
    new_eq = eq_lookup(new_inchikey)

    def fmt_eq(e):
        return f"{e[1]:.2f} via {e[0]}" if e else "(no MetaNetX/eQ hit for deprotonated InChIKey)"

    eq_changed = (cur_eq != new_eq)
    print(f"\nΔG (eQuilibrator):")
    print(f"  current : {fmt_eq(cur_eq)}")
    print(f"  new     : {fmt_eq(new_eq)}{' [CHANGE]' if eq_changed else ''}")

    # ---- Group Contribution (not recomputed) ---------------------------
    cur_gc = cur_thermo.get('Group contribution') if isinstance(cur_thermo, dict) else None
    print(f"\nΔG (Group Contribution):")
    print(f"  current : {cur_gc}")
    print(f"  new     : (would require MFAToolkit re-run on the new mol file — not computed here)")

    # ---- pKa attribution (informational) -------------------------------
    aliases = helper.loadMSAliases()
    cpd_aliases = aliases.get(cpd_id, {})
    pka_chain = []
    for db in ['KEGG', 'MetaCyc']:
        for ext in cpd_aliases.get(db, []):
            pka_chain.append(f"{db}:{ext}")
    print(f"\npKa / pKb (alias-based, NOT structure-derived):")
    print(f"  current : pka={cpd.get('pka')}  pkb={cpd.get('pkb')}")
    if pka_chain:
        print(f"  chain   : {' → '.join(pka_chain[:6])}{'…' if len(pka_chain) > 6 else ''}")
    print(f"  note    : updating the structure does NOT change pka unless the alias")
    print(f"            mapping also changes.")

    # ---- Summary -------------------------------------------------------
    changed = []
    if new_formula != cpd.get('formula'):    changed.append('formula')
    if str(new_charge) != str(cpd.get('charge')):  changed.append('charge')
    if new_inchikey != cpd.get('inchikey'):  changed.append('inchikey')
    if eq_changed:                           changed.append('deltag (eQ)')
    print(f"\nSummary: {len(changed)} field(s) would change: {', '.join(changed) if changed else '(none)'}")
    print(f"         GC ΔG and pka/pkb attribution unaffected by this preview "
          f"(re-run the GC + pKa update steps after committing the new structure).")


if __name__ == '__main__':
    main()
