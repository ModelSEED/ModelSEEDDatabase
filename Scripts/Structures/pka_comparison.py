"""Compare protonation corrections against ChemAxon pKa values.

For each compound where a protonation correction was made (or a
protonation difference with PubChem was detected), this script checks
whether the ChemAxon pKa values from the ModelSEED database support
the stored charge, the OpenBabel-corrected charge, or the PubChem
charge.

Logic:
  - pKa < 7  →  group is deprotonated at pH 7  (loses H, charge −1)
  - pKa > 7  →  group stays protonated at pH 7  (keeps H, neutral)
  - pKb > 7  →  basic group is protonated at pH 7 (gains H, charge +1)
  - pKb < 7  →  basic group is neutral at pH 7
"""

import csv
import glob
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DB_DIR = os.path.join(os.path.dirname(SCRIPT_DIR),
                      "ModelSEEDDatabase", "Biochemistry")


def load_pka_from_db():
    """Load pKa and pKb values from compound_XX.tsv files.

    Returns dict: {cpd_id: {'pka': [(frag, atom, value), ...],
                             'pkb': [(frag, atom, value), ...],
                             'name': str, 'charge': int}}
    """
    pka_data = {}
    compound_files = sorted(glob.glob(os.path.join(DB_DIR, "compound_*.tsv")))
    for filepath in compound_files:
        with open(filepath) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                cpd_id = row['id']
                pka_str = row.get('pka', '') or ''
                pkb_str = row.get('pkb', '') or ''
                name = row.get('name', '')
                charge = row.get('charge', '0')

                pka_vals = _parse_pka_string(pka_str)
                pkb_vals = _parse_pka_string(pkb_str)

                if pka_vals or pkb_vals:
                    pka_data[cpd_id] = {
                        'pka': pka_vals,
                        'pkb': pkb_vals,
                        'name': name,
                        'db_charge': int(charge) if charge else 0,
                    }
    return pka_data


def _parse_pka_string(s):
    """Parse '1:14:12.60;1:22:3.33' → [(1, 14, 12.60), (1, 22, 3.33)]."""
    if not s or s == 'null':
        return []
    entries = []
    for part in s.split(';'):
        fields = part.split(':')
        if len(fields) == 3:
            try:
                entries.append((int(fields[0]), int(fields[1]),
                                float(fields[2])))
            except ValueError:
                continue
    return entries


def predict_charge_from_pka(pka_vals, pkb_vals, ph=7.0):
    """Predict net charge contribution from ionizable groups at given pH.

    Returns:
      deprotonations: number of acidic groups deprotonated (pKa < pH)
      protonations: number of basic groups protonated (pKb > pH)
      net_charge_contribution: -deprotonations + protonations
    """
    deprotonations = sum(1 for _, _, v in pka_vals if v < ph)
    protonations = sum(1 for _, _, v in pkb_vals if v > ph)
    return deprotonations, protonations, -deprotonations + protonations


def main():
    print("=" * 80)
    print("pKa-Based Protonation Analysis")
    print("Comparing stored/corrected charges against ChemAxon pKa predictions")
    print("=" * 80)

    pka_data = load_pka_from_db()
    print(f"\nLoaded pKa data for {len(pka_data)} compounds from database.\n")

    # --- Section 1: Compounds with pH 7 corrections (Phase 4) ---
    corrections_file = os.path.join(SCRIPT_DIR,
                                    "pubchem_protonation_corrections.tsv")
    corrections = []
    if os.path.exists(corrections_file):
        with open(corrections_file) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            corrections = list(reader)

    print("-" * 80)
    print("SECTION 1: OpenBabel pH 7 Corrections vs ChemAxon pKa")
    print("-" * 80)

    if not corrections:
        print("  No protonation corrections found.\n")
    else:
        for row in corrections:
            cpd_id = row['cpd_id']
            name = row['compound_name']
            stored_charge = int(row['stored_charge'])
            ph7_charge = int(row['ph7_charge'])

            print(f"\n  {cpd_id} ({name})")
            print(f"    Stored charge:            {stored_charge}")
            print(f"    OpenBabel pH 7 charge:    {ph7_charge}")

            info = pka_data.get(cpd_id)
            if not info:
                print("    ChemAxon pKa:             NOT AVAILABLE")
                continue

            pka_vals = info['pka']
            pkb_vals = info['pkb']

            # Display all pKa values
            pka_sorted = sorted(pka_vals, key=lambda x: x[2])
            print(f"    pKa values:               "
                  f"{', '.join(f'{v:.2f} (atom {a})' for _, a, v in pka_sorted)}")
            if pkb_vals:
                pkb_sorted = sorted(pkb_vals, key=lambda x: x[2])
                print(f"    pKb values:               "
                      f"{', '.join(f'{v:.2f} (atom {a})' for _, a, v in pkb_sorted)}")

            # Classify each group
            deprot, prot, net = predict_charge_from_pka(pka_vals, pkb_vals)
            print(f"    Acidic groups w/ pKa < 7: {deprot} "
                  f"(deprotonated → charge −{deprot})")
            print(f"    Basic groups w/ pKb > 7:  {prot} "
                  f"(protonated → charge +{prot})")
            print(f"    pKa-predicted net charge:  {net}")

            # Verdict
            if net == stored_charge:
                print(f"    >>> pKa SUPPORTS stored charge ({stored_charge}), "
                      f"NOT the OpenBabel correction ({ph7_charge})")
            elif net == ph7_charge:
                print(f"    >>> pKa SUPPORTS OpenBabel correction ({ph7_charge}), "
                      f"NOT stored charge ({stored_charge})")
            else:
                print(f"    >>> pKa predicts {net}, which matches NEITHER "
                      f"stored ({stored_charge}) nor OpenBabel ({ph7_charge})")

            # Flag borderline pKa values (6.0-8.0)
            borderline = [(a, v) for _, a, v in pka_vals if 6.0 <= v <= 8.0]
            if borderline:
                print(f"    *** BORDERLINE pKa: "
                      f"{', '.join(f'atom {a} pKa={v:.2f}' for a, v in borderline)} "
                      f"— near pH 7, result is method-dependent")

    # --- Section 2: PubChem protonation differences ---
    diffs_file = os.path.join(SCRIPT_DIR, "pubchem_protonation_diffs.tsv")
    diffs = []
    if os.path.exists(diffs_file):
        with open(diffs_file) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            diffs = list(reader)

    print("\n" + "-" * 80)
    print("SECTION 2: PubChem Protonation Differences vs ChemAxon pKa")
    print("-" * 80)

    if not diffs:
        print("  No protonation differences found.\n")
    else:
        for row in diffs:
            cpd_id = row['cpd_id']
            name = row['compound_name']
            stored_ik = row['stored_inchikey']
            pubchem_ik = row['pubchem_inchikey']

            # Determine charges from SMILES
            from rdkit import Chem
            stored_mol = Chem.MolFromSmiles(row['stored_smiles'])
            pubchem_mol = Chem.MolFromSmiles(row['pubchem_smiles'])
            stored_q = Chem.GetFormalCharge(stored_mol) if stored_mol else '?'
            pubchem_q = Chem.GetFormalCharge(pubchem_mol) if pubchem_mol else '?'

            print(f"\n  {cpd_id} ({name})")
            print(f"    Stored InChIKey suffix:    ...{stored_ik[-1]}"
                  f"  (charge: {stored_q})")
            print(f"    PubChem InChIKey suffix:   ...{pubchem_ik[-1]}"
                  f"  (charge: {pubchem_q})")

            info = pka_data.get(cpd_id)
            if not info:
                print("    ChemAxon pKa:             NOT AVAILABLE")
                continue

            pka_vals = info['pka']
            pkb_vals = info['pkb']

            pka_sorted = sorted(pka_vals, key=lambda x: x[2])
            print(f"    pKa values:               "
                  f"{', '.join(f'{v:.2f} (atom {a})' for _, a, v in pka_sorted)}")
            if pkb_vals:
                pkb_sorted = sorted(pkb_vals, key=lambda x: x[2])
                print(f"    pKb values:               "
                      f"{', '.join(f'{v:.2f} (atom {a})' for _, a, v in pkb_sorted)}")

            deprot, prot, net = predict_charge_from_pka(pka_vals, pkb_vals)
            print(f"    Acidic groups w/ pKa < 7: {deprot}")
            print(f"    Basic groups w/ pKb > 7:  {prot}")
            print(f"    pKa-predicted net charge:  {net}")

            if net == stored_q:
                print(f"    >>> pKa SUPPORTS stored charge ({stored_q}), "
                      f"NOT PubChem ({pubchem_q})")
            elif net == pubchem_q:
                print(f"    >>> pKa SUPPORTS PubChem charge ({pubchem_q}), "
                      f"NOT stored ({stored_q})")
            else:
                print(f"    >>> pKa predicts {net}, which matches NEITHER "
                      f"stored ({stored_q}) nor PubChem ({pubchem_q})")

            borderline = [(a, v) for _, a, v in pka_vals if 6.0 <= v <= 8.0]
            if borderline:
                print(f"    *** BORDERLINE pKa: "
                      f"{', '.join(f'atom {a} pKa={v:.2f}' for a, v in borderline)} "
                      f"— near pH 7, result is method-dependent")

    # --- Section 3: Summary ---
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Broad stats: how many compounds have pKa data, how many have
    # borderline values, how many have all groups clearly below/above 7
    total_with_pka = len(pka_data)
    borderline_compounds = 0
    for cpd_id, info in pka_data.items():
        for _, _, v in info['pka']:
            if 6.0 <= v <= 8.0:
                borderline_compounds += 1
                break
    print(f"\n  Compounds with pKa data:    {total_with_pka}")
    print(f"  Compounds with borderline"
          f"\n    pKa (6.0-8.0):            {borderline_compounds}")

    # For the corrected compounds, tally verdicts
    supported_stored = 0
    supported_corrected = 0
    ambiguous = 0
    for row in corrections:
        cpd_id = row['cpd_id']
        info = pka_data.get(cpd_id)
        if not info:
            continue
        stored_charge = int(row['stored_charge'])
        ph7_charge = int(row['ph7_charge'])
        _, _, net = predict_charge_from_pka(info['pka'], info['pkb'])
        if net == stored_charge:
            supported_stored += 1
        elif net == ph7_charge:
            supported_corrected += 1
        else:
            ambiguous += 1

    print(f"\n  pH 7 corrections (Phase 4):")
    print(f"    pKa supports STORED charge:     {supported_stored}")
    print(f"    pKa supports CORRECTED charge:  {supported_corrected}")
    print(f"    pKa supports NEITHER:           {ambiguous}")


if __name__ == '__main__':
    main()
