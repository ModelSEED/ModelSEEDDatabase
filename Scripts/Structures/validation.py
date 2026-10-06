import sys
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem.inchi import MolFromInchi, MolToInchi, InchiToInchiKey


def validate_compound(args):
    """Validate a single compound's SMILES, InChI, and InChIKey consistency.

    Returns (cpd_id, issues_list).
    """
    cpd_id, smiles, inchikey, inchi = args
    issues = []

    # Parse SMILES
    smi_mol = None
    if smiles:
        smi_mol = Chem.MolFromSmiles(smiles)
        if smi_mol is None:
            issues.append(f"SMILES fails RDKit parsing: {smiles[:80]}")

    # Parse InChI
    inchi_mol = None
    if inchi:
        inchi_mol = MolFromInchi(inchi)
        if inchi_mol is None:
            issues.append(f"InChI fails RDKit parsing: {inchi[:80]}")

    # Recompute InChIKey from InChI and compare
    if inchi_mol and inchikey:
        try:
            computed_key = InchiToInchiKey(inchi)
            if computed_key != inchikey:
                issues.append(
                    f"InChIKey mismatch: stored={inchikey} computed={computed_key}"
                )
        except Exception as e:
            issues.append(f"InChIKey computation failed: {e}")

    # Cross-check SMILES vs InChI round-trip
    if smi_mol and inchi_mol:
        try:
            inchi_from_smi = MolToInchi(smi_mol)
            key_from_smi = InchiToInchiKey(inchi_from_smi) if inchi_from_smi else None
            key_from_inchi = InchiToInchiKey(inchi) if inchi else None
            if key_from_smi and key_from_inchi and key_from_smi != key_from_inchi:
                issues.append(
                    f"SMILES/InChI disagree: key_from_SMILES={key_from_smi} key_from_InChI={key_from_inchi}"
                )
        except Exception as e:
            issues.append(f"Cross-validation failed: {e}")

    # Flag compounds missing all structure types
    if not smiles and not inchi and not inchikey:
        issues.append("No structure data at all")

    return cpd_id, issues


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else 'Unique_ModelSEED_Structures_new.txt'
    num_workers = min(cpu_count(), 64)

    print(f"Loading {input_path}...")

    # Group rows by compound ID
    compounds = defaultdict(lambda: {'SMILE': None, 'InChIKey': None, 'InChI': None})
    row_count = 0

    with open(input_path, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 6:
                continue
            row_count += 1
            cpd_id, typ, aliases, formula, charge, struct = parts
            if typ in ('SMILE', 'InChIKey', 'InChI'):
                compounds[cpd_id][typ] = struct

    print(f"  Rows: {row_count}")
    print(f"  Compounds: {len(compounds)}")

    # Prepare work items
    work_items = [
        (cpd_id, d['SMILE'], d['InChIKey'], d['InChI'])
        for cpd_id, d in sorted(compounds.items())
    ]

    # Validate in parallel
    print(f"\nValidating with {num_workers} workers...")
    total_issues = 0
    failed_compounds = []

    with Pool(num_workers) as pool:
        for cpd_id, issues in pool.imap_unordered(validate_compound, work_items, chunksize=256):
            if issues:
                total_issues += len(issues)
                failed_compounds.append((cpd_id, issues))

    # Report
    failed_compounds.sort()
    print(f"\nResults:")
    print(f"  Compounds checked: {len(compounds)}")
    print(f"  Compounds with issues: {len(failed_compounds)}")
    print(f"  Total issues: {total_issues}")

    if failed_compounds:
        # Categorize issues
        stereo_mismatch = []   # Same connectivity, different stereochemistry
        diff_compound = []     # Completely different first block of InChIKey
        inchi_parse_fail = []  # InChI fails RDKit parsing
        smiles_parse_fail = [] # SMILES fails RDKit parsing
        other = []

        for cpd_id, issues in failed_compounds:
            for issue in issues:
                if "InChI fails RDKit parsing" in issue:
                    inchi_parse_fail.append((cpd_id, issue))
                elif "SMILES fails RDKit parsing" in issue:
                    smiles_parse_fail.append((cpd_id, issue))
                elif "SMILES/InChI disagree" in issue:
                    # Extract keys and compare first block (connectivity)
                    parts = issue.split('=')
                    key_smi = parts[1].split()[0]
                    key_inchi = parts[2].strip()
                    if key_smi.split('-')[0] == key_inchi.split('-')[0]:
                        stereo_mismatch.append((cpd_id, issue))
                    else:
                        diff_compound.append((cpd_id, issue))
                elif "InChIKey mismatch" in issue:
                    stereo_mismatch.append((cpd_id, issue))
                else:
                    other.append((cpd_id, issue))

        print(f"\n--- Issue Categories ---")
        print(f"\n[1] Stereochemistry mismatch (same connectivity, different stereo): {len(stereo_mismatch)}")
        for cpd_id, issue in stereo_mismatch:
            print(f"  {cpd_id}: {issue}")

        print(f"\n[2] Different compounds (SMILES and InChI encode different molecules): {len(diff_compound)}")
        for cpd_id, issue in diff_compound:
            print(f"  {cpd_id}: {issue}")

        print(f"\n[3] InChI fails RDKit parsing: {len(inchi_parse_fail)}")
        for cpd_id, issue in inchi_parse_fail:
            print(f"  {cpd_id}: {issue}")

        print(f"\n[4] SMILES fails RDKit parsing: {len(smiles_parse_fail)}")
        for cpd_id, issue in smiles_parse_fail:
            print(f"  {cpd_id}: {issue}")

        if other:
            print(f"\n[5] Other issues: {len(other)}")
            for cpd_id, issue in other:
                print(f"  {cpd_id}: {issue}")
    else:
        print("\nAll compounds passed validation.")


if __name__ == '__main__':
    main()
