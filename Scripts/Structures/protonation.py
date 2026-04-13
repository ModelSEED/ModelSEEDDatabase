"""pH 7 protonation normalization for ModelSEED compounds.

Uses OpenBabel's CorrectForPH(7.0) to normalize protonation states.
Checks all compounds and corrects only those not already at pH 7.
"""

import csv
import logging
import os
import sys
from multiprocessing import Pool, cpu_count

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(BASE_DIR, "ModelSEEDDatabase", "Libs", "Python"))

from rdkit import Chem

from structure_compare import compare_inchikeys, compute_formula_charge_from_inchi

logger = logging.getLogger("pubchem_validate")


def _is_over_deprotonated(stored_smiles, ph7_smiles):
    """Check if OpenBabel over-deprotonated a free phosphate group.

    OpenBabel's CorrectForPH(7.0) treats all phosphate -OH groups the
    same, but free (non-ester) polyphosphates like pyrophosphate have
    pKa4 ~9.4 — well above pH 7.  In contrast, ester-linked phosphates
    (nucleotide di/triphosphates) have terminal pKa ~6.5 and ARE
    correctly deprotonated at pH 7.

    This function checks whether ALL phosphorus atoms in the molecule
    are free (no C-O-P ester bond on any P in the chain).  Only in that
    case can the last terminal -OH have pKa ~9.3.

    Returns True if the correction should be rejected.
    """
    stored_mol = Chem.MolFromSmiles(stored_smiles) if stored_smiles else None
    ph7_mol = Chem.MolFromSmiles(ph7_smiles) if ph7_smiles else None
    if stored_mol is None or ph7_mol is None:
        return False

    stored_chg = Chem.GetFormalCharge(stored_mol)
    ph7_chg = Chem.GetFormalCharge(ph7_mol)

    # Only check cases where charge became more negative
    if ph7_chg >= stored_chg:
        return False

    # Collect all P atoms and check if ANY has a C-O-P ester bond.
    # If even one P in the molecule is ester-linked, the polyphosphate
    # chain is attached to an organic moiety and all terminal pKa values
    # are ~6.5 (correctly deprotonated at pH 7).
    p_atoms = [a for a in ph7_mol.GetAtoms() if a.GetAtomicNum() == 15]
    if not p_atoms:
        return False

    chain_has_ester = False
    for p_atom in p_atoms:
        for neighbor in p_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # oxygen
                for nn in neighbor.GetNeighbors():
                    if nn.GetAtomicNum() == 6:  # carbon — ester bond
                        chain_has_ester = True
                        break
            if chain_has_ester:
                break
        if chain_has_ester:
            break

    if chain_has_ester:
        # Ester-linked phosphate chain (e.g. ATP, ADP): terminal pKa
        # ~6.5, deprotonated at pH 7.  OpenBabel is correct.
        return False

    # All P atoms are free (no ester).  Check if any terminal P has
    # both non-bridging oxygens fully deprotonated — that means the
    # pKa4 (~9.3) proton was removed, which is wrong at pH 7.
    for p_atom in p_atoms:
        neg_oxygens = 0
        for neighbor in p_atom.GetNeighbors():
            if (neighbor.GetAtomicNum() == 8
                    and neighbor.GetFormalCharge() == -1
                    and neighbor.GetDegree() == 1):
                neg_oxygens += 1
        if neg_oxygens >= 2:
            return True

    return False


# Module-level handles for OpenBabel (set once per worker process)
_pybel = None
_openbabel = None


def _ensure_openbabel():
    """Initialize OpenBabel if not already done (for non-pool usage)."""
    global _pybel, _openbabel
    if _pybel is not None:
        return
    from openbabel import openbabel as _ob
    _ob.obErrorLog.SetOutputLevel(_ob.obError)
    _ob.obErrorLog.StopLogging()
    from openbabel import pybel as _pb
    _pybel = _pb
    _openbabel = _ob


def _init_openbabel_worker():
    """Pool initializer: import OpenBabel once per worker process."""
    _ensure_openbabel()


def normalize_to_ph7(smiles):
    """Normalize a SMILES string to pH 7 protonation using OpenBabel.

    Returns dict with {smiles, inchi, inchikey} or None on failure.
    When called inside a multiprocessing Pool initialized with
    _init_openbabel_worker, uses pre-imported handles. Otherwise
    initializes OpenBabel on first call.
    """
    if not smiles or not smiles.strip():
        return None
    _ensure_openbabel()
    try:
        mol = _pybel.readstring('smi', smiles)
        mol.OBMol.CorrectForPH(7.0)

        new_smiles = mol.write('smi').strip().split()[0]
        new_inchi = mol.write('inchi').strip()
        new_inchikey = mol.write('inchikey').strip()

        if not new_inchi or not new_inchikey:
            return None

        return {
            'smiles': new_smiles,
            'inchi': new_inchi,
            'inchikey': new_inchikey,
        }
    except Exception:
        return None


def _ph7_normalize_worker(args):
    """Worker for parallel pH 7 normalization (module-level for pickling)."""
    cpd_id, smiles = args
    result = normalize_to_ph7(smiles)
    return cpd_id, result


def run_phase4_protonation(candidates, structures, names=None,
                           corrections_file=None):
    """Phase 4: Normalize all compounds to pH 7 protonation.

    For each compound with a stored SMILES:
    1. Compute pH 7 form via OpenBabel CorrectForPH(7.0)
    2. Compare pH 7 InChIKey to stored InChIKey
    3. If they match -> skip (already at pH 7)
    4. If protonation differs -> correct to pH 7 form
    5. If connectivity differs -> skip (OpenBabel artifact)

    Returns dict of corrections compatible with apply_corrections():
        {cpd_id: {smiles, inchi, inchikey, result_type, strategy, query}}
    """
    if corrections_file is None:
        corrections_file = os.path.join(SCRIPT_DIR,
                                        "pubchem_protonation_corrections.tsv")
    # CPU-bound work: use cpu_count, not the network-worker count
    num_cpus = min(cpu_count(), 32)

    # Collect all compounds with SMILES
    work_items = []
    for cpd_id in candidates:
        stored = structures.get(cpd_id, {})
        smiles = stored.get("smiles", "")
        if smiles:
            work_items.append((cpd_id, smiles))

    logger.info("Phase 4: pH 7 protonation normalization")
    logger.info("  Compounds with SMILES: %d", len(work_items))

    if not work_items:
        logger.info("  No compounds to normalize.")
        return {}

    # Normalize in parallel
    ph7_results = {}  # cpd_id -> {smiles, inchi, inchikey} or None

    try:
        from tqdm import tqdm as _tqdm
    except ImportError:
        _tqdm = None

    if _tqdm:
        pbar = _tqdm(total=len(work_items), desc="Phase 4: pH 7 normalization")

    done = 0
    with Pool(num_cpus, initializer=_init_openbabel_worker) as pool:
        for cpd_id, result in pool.imap_unordered(
                _ph7_normalize_worker, work_items, chunksize=64):
            ph7_results[cpd_id] = result
            done += 1
            if _tqdm:
                pbar.update(1)
            elif done % 5000 == 0:
                logger.info("  Phase 4 progress: %d/%d", done, len(work_items))

    if _tqdm:
        pbar.close()

    # Classify results
    corrections = {}
    stats = {
        'already_correct': 0,
        'protonation_corrected': 0,
        'connectivity_changed': 0,
        'normalization_failed': 0,
    }
    correction_rows = []  # for TSV output

    for cpd_id in candidates:
        stored = structures.get(cpd_id, {})
        stored_ik = stored.get("inchikey", "")
        stored_smi = stored.get("smiles", "")
        if not stored_smi or cpd_id not in ph7_results:
            continue

        ph7 = ph7_results[cpd_id]
        if ph7 is None:
            stats['normalization_failed'] += 1
            continue

        ph7_ik = ph7['inchikey']
        result_type = compare_inchikeys(stored_ik, ph7_ik)

        if result_type == "MATCH":
            stats['already_correct'] += 1
            continue

        if result_type == "PROTONATION_DIFF":
            # Reject if OpenBabel over-deprotonated a free phosphate
            if _is_over_deprotonated(stored_smi, ph7['smiles']):
                stats['connectivity_changed'] += 1
                logger.info("  %s: pH 7 correction rejected "
                            "(over-deprotonated free phosphate, pKa > 7)",
                            cpd_id)
                continue
            # pH 7 form differs only in protonation -> apply correction
            formula, charge = compute_formula_charge_from_inchi(ph7['inchi'])
            corrections[cpd_id] = {
                'smiles': ph7['smiles'],
                'inchi': ph7['inchi'],
                'inchikey': ph7['inchikey'],
                'pubchem_cid': 'pH7',
                'result_type': 'PH7_CORRECTION',
                'strategy': 'ph7_normalization',
                'query': 'OpenBabel_CorrectForPH(7.0)',
                'validation_reason': (
                    f"Protonation normalized to pH 7 "
                    f"(stored={stored_ik[-1]}, ph7={ph7_ik[-1]})"),
            }
            cpd_names = names.get(cpd_id, []) if names else []
            name_str = cpd_names[0] if cpd_names else ""
            correction_rows.append({
                'cpd_id': cpd_id,
                'compound_name': name_str,
                'action': 'corrected',
                'stored_inchikey': stored_ik,
                'ph7_inchikey': ph7_ik,
                'stored_smiles': stored_smi,
                'ph7_smiles': ph7['smiles'],
                'stored_formula': stored.get('formula', ''),
                'ph7_formula': formula or '',
                'stored_charge': stored.get('charge', ''),
                'ph7_charge': str(charge) if charge is not None else '',
            })
            stats['protonation_corrected'] += 1

        elif result_type in ("STEREO_DIFF", "MISMATCH"):
            # Connectivity or stereo changed — OpenBabel artifact, skip
            stats['connectivity_changed'] += 1

        else:
            stats['normalization_failed'] += 1

    # Write corrections TSV
    if correction_rows:
        fieldnames = ['cpd_id', 'compound_name', 'action',
                      'stored_inchikey', 'ph7_inchikey',
                      'stored_smiles', 'ph7_smiles',
                      'stored_formula', 'ph7_formula',
                      'stored_charge', 'ph7_charge']
        with open(corrections_file, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in sorted(correction_rows, key=lambda r: r['cpd_id']):
                writer.writerow(row)
        logger.info("  Protonation corrections log: %s", corrections_file)

    logger.info("  Already at pH 7: %d", stats['already_correct'])
    logger.info("  Protonation corrected: %d", stats['protonation_corrected'])
    logger.info("  Connectivity changed (skipped): %d",
                stats['connectivity_changed'])
    logger.info("  Normalization failed: %d", stats['normalization_failed'])

    return corrections
