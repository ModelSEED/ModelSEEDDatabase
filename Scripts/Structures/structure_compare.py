"""Structure comparison and classification utilities for PubChem validation.

Provides functions for comparing InChIKeys, classifying mismatches,
validating stereo corrections, and related structural analysis.
"""

import logging
import os
import re
import sys

from rdkit import Chem, DataStructs
from rdkit.Chem import FindPotentialStereo, StereoSpecified
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(BASE_DIR, "ModelSEEDDatabase", "Libs", "Python"))
from BiochemPy import Compounds, InChIs

TANIMOTO_CUTOFF = 0.8
MAX_NAMES = 3

_taut_enum = TautomerEnumerator()
logger = logging.getLogger("pubchem_validate")


def pick_best_names(name_list, max_names=MAX_NAMES):
    """Pick the best names for PubChem lookup.

    Returns up to max_names names sorted by quality score (best first).
    Prefers longer, descriptive names starting with letters.
    """
    cleaned = []
    for n in name_list:
        n = re.sub(r"<[^>]+>", "", n)
        n = n.strip()
        if n:
            cleaned.append(n)
    if not cleaned:
        return []
    scored = []
    for n in cleaned:
        score = len(n)
        if re.match(r"^[A-Za-z]", n):
            score += 100
        if re.search(r"[a-z]", n):
            score += 50
        if re.match(r"^[A-Z0-9+\-\s]{1,5}$", n):
            score -= 200
        scored.append((score, n))
    scored.sort(reverse=True)
    seen = set()
    result = []
    for _, n in scored:
        key = n.lower()
        if key not in seen:
            seen.add(key)
            result.append(n)
            if len(result) >= max_names:
                break
    return result


# ---------------------------------------------------------------------------
# InChIKey comparison
# ---------------------------------------------------------------------------

def compare_inchikeys(stored, pubchem):
    """Compare two InChIKeys and return classification."""
    if not stored or not pubchem:
        return "NO_KEY"
    if stored == pubchem:
        return "MATCH"
    s_parts = stored.split("-")
    p_parts = pubchem.split("-")
    if len(s_parts) != 3 or len(p_parts) != 3:
        return "MISMATCH"
    if s_parts[0] == p_parts[0] and s_parts[1] == p_parts[1]:
        return "PROTONATION_DIFF"
    if s_parts[0] == p_parts[0]:
        return "STEREO_DIFF"
    return "MISMATCH"


def compare_inchi_layers(stored_inchi, pubchem_inchi):
    """Compare InChI strings at different layer levels for finer classification.

    Returns dict with boolean flags for connectivity, stereo, and protonation
    layer matches, or None if either InChI is missing/unparseable.
    """
    if not stored_inchi or not pubchem_inchi:
        return None
    try:
        s_formula, s_layers = InChIs.parse(stored_inchi)
        p_formula, p_layers = InChIs.parse(pubchem_inchi)
    except Exception:
        return None

    # Compare without stereo layers
    s_no_stereo = InChIs.build(s_formula, s_layers,
                               remove=('b', 't', 'm', 's'))
    p_no_stereo = InChIs.build(p_formula, p_layers,
                               remove=('b', 't', 'm', 's'))

    # Compare without protonation
    s_no_prot = InChIs.build(s_formula, s_layers, remove=('p', 'q'))
    p_no_prot = InChIs.build(p_formula, p_layers, remove=('p', 'q'))

    # Compare connectivity only (no stereo, no protonation)
    s_conn = InChIs.build(s_formula, s_layers,
                          remove=('b', 't', 'm', 's', 'p', 'q'))
    p_conn = InChIs.build(p_formula, p_layers,
                          remove=('b', 't', 'm', 's', 'p', 'q'))

    return {
        "connectivity_match": s_conn == p_conn and s_conn != "",
        "stereo_match": s_no_prot == p_no_prot and s_no_prot != "",
        "protonation_match": s_no_stereo == p_no_stereo and s_no_stereo != "",
    }


# ---------------------------------------------------------------------------
# Correction validation helpers
# ---------------------------------------------------------------------------

def count_defined_stereo(smiles):
    """Count defined stereocenters and E/Z bonds in a SMILES string.

    Returns (specified_count, potential_count) or (None, None) on failure.
    """
    if not smiles:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    stereo_info = FindPotentialStereo(mol)
    potential = len(stereo_info)
    specified = sum(1 for s in stereo_info
                    if s.specified == StereoSpecified.Specified)
    return specified, potential


def _check_stereo_compatibility(stored_struct, pubchem_struct):
    """Verify that stereocenters defined in both molecules agree.

    Uses InChI stereo layers (t-layer for tetrahedral, b-layer for double
    bond) which use canonical atom numbering, avoiding the pitfalls of
    SMILES-based atom mapping.

    For STEREO_DIFF corrections: if PubChem inverts existing stereocenters,
    the correction should be rejected.

    Returns (compatible: bool, detail: str).
    """
    stored_inchi = stored_struct.get("inchi", "")
    pub_inchi = pubchem_struct.get("inchi", "")
    if not stored_inchi or not pub_inchi:
        return True, "skipped (no InChI data)"

    try:
        _, s_layers = InChIs.parse(stored_inchi)
        _, p_layers = InChIs.parse(pub_inchi)
    except Exception:
        return True, "skipped (InChI parse error)"

    # Compare tetrahedral stereo (t-layer)
    s_t = s_layers.get('t', '')
    p_t = p_layers.get('t', '')

    # Compare double bond stereo (b-layer)
    s_b = s_layers.get('b', '')
    p_b = p_layers.get('b', '')

    if not s_t and not s_b:
        # Stored has no stereo defined in InChI — any PubChem stereo is new
        return True, "compatible (stored has no InChI stereo layers)"

    inversions = 0
    checked = 0

    # Parse t-layer: format like "+1,+3,-5,?" where numbers are atom indices
    # and +/- indicates configuration. '?' means undefined.
    if s_t and p_t:
        s_centers = {}  # atom_num -> config (+/-)
        for part in re.findall(r'([+-]?\d+)', s_t):
            if part.startswith('+') or part.startswith('-'):
                s_centers[part[1:]] = part[0]
        p_centers = {}
        for part in re.findall(r'([+-]?\d+)', p_t):
            if part.startswith('+') or part.startswith('-'):
                p_centers[part[1:]] = part[0]
        # Check shared defined centers
        for atom_num, s_config in s_centers.items():
            if atom_num in p_centers:
                checked += 1
                if s_config != p_centers[atom_num]:
                    inversions += 1

    # Parse b-layer: format like "1-3+,5-7-" for double bond stereo
    if s_b and p_b:
        s_bonds = {}  # bond_spec -> config (+/-)
        for match in re.finditer(r'(\d+-\d+)([+-])', s_b):
            s_bonds[match.group(1)] = match.group(2)
        p_bonds = {}
        for match in re.finditer(r'(\d+-\d+)([+-])', p_b):
            p_bonds[match.group(1)] = match.group(2)
        for bond_spec, s_config in s_bonds.items():
            if bond_spec in p_bonds:
                checked += 1
                if s_config != p_bonds[bond_spec]:
                    inversions += 1

    if inversions > 0:
        return False, (f"stereo_inversion: {inversions} of {checked} "
                       f"shared stereocenters have different configuration")
    return True, f"compatible ({checked} shared stereocenters agree)"


def compute_formula_charge_from_inchi(inchi_str):
    """Compute molecular formula and charge from an InChI string.

    Uses BiochemPy InChIs utilities to parse layers and adjust for protons.
    Returns (formula, charge) or (None, None) on failure.
    """
    if not inchi_str:
        return None, None
    try:
        formula, layers = InChIs.parse(inchi_str)
        chg = InChIs.charge(layers['q'], layers['p'])
        # Merge multi-component formulas (e.g. "C10H12N5O3.Co") before
        # adjusting protons, since adjust_protons requires a single formula
        if '.' in formula:
            formula, _ = Compounds.mergeFormula(formula)
        if layers['p']:
            formula, _ = InChIs.adjust_protons(formula, layers['p'])
        return formula, chg
    except Exception:
        return None, None


def _verify_smiles_inchikey_consistency(smiles, expected_inchikey):
    """Check that SMILES generates the same InChIKey as PubChem reports.

    Returns (ok: bool, detail: str). Compares only the first (connectivity)
    block of the InChIKey, since protonation/stereo differences between
    SMILES and InChI representations are expected in some cases.
    """
    if not smiles or not expected_inchikey:
        return True, "skipped (missing data)"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "cannot parse SMILES"
    derived_inchikey = Chem.inchi.MolToInchiKey(mol)
    if not derived_inchikey:
        return False, "cannot generate InChIKey from SMILES"
    # Compare first block (connectivity) — stereo/protonation may differ
    # due to RDKit vs PubChem normalization differences
    d1 = derived_inchikey.split("-")[0]
    e1 = expected_inchikey.split("-")[0]
    if d1 != e1:
        return False, (f"InChIKey mismatch: SMILES->{derived_inchikey[:14]} "
                       f"vs reported {expected_inchikey[:14]}")
    return True, "consistent"


def _classify_mismatch(stored_struct, pubchem_struct):
    """Classify a MISMATCH into a sub-category for the review report.

    Returns a diagnostic string describing why the structures differ.
    Used for manual review — MISMATCH is never auto-corrected.
    """
    stored_smiles = stored_struct.get("smiles", "")
    pub_smiles = pubchem_struct.get("smiles", "")
    stored_formula = stored_struct.get("formula", "")
    pub_inchi = pubchem_struct.get("inchi", "")

    pub_formula, _ = compute_formula_charge_from_inchi(pub_inchi)
    # Fall back to formula from pub_struct if InChI-derived formula unavailable
    if not pub_formula:
        pub_formula = pubchem_struct.get("formula", "")
    sf = stored_formula.strip() if stored_formula else ""
    pf = pub_formula.strip() if pub_formula else ""

    if sf and pf and sf != pf:
        # Normalize dot-separated salt formulas (e.g. "6CN.Fe" -> "C6FeN6")
        sf_merged, _ = Compounds.mergeFormula(sf)
        pf_merged, _ = Compounds.mergeFormula(pf)
        if sf_merged != "null" and pf_merged != "null" and sf_merged != pf_merged:
            # Check if difference is only in hydrogen count (protonation)
            sf_atoms = Compounds.parseFormula(sf_merged)
            pf_atoms = Compounds.parseFormula(pf_merged)
            sf_no_h = {k: v for k, v in sf_atoms.items() if k != "H"}
            pf_no_h = {k: v for k, v in pf_atoms.items() if k != "H"}
            if sf_no_h == pf_no_h:
                h_diff = pf_atoms.get("H", 0) - sf_atoms.get("H", 0)
                return (f"IGNORE:protonation_formula_diff (formula: "
                        f"stored={sf}, pubchem={pf}, H_diff={h_diff:+d})")
            return (f"IGNORE:wrong_mapping (formula: stored={sf}, "
                    f"pubchem={pf})")

    # InChI layer comparison for finer classification
    stored_inchi = stored_struct.get("inchi", "")
    if stored_inchi and pub_inchi:
        layer_cmp = compare_inchi_layers(stored_inchi, pub_inchi)
        if layer_cmp and layer_cmp["connectivity_match"]:
            if layer_cmp["stereo_match"]:
                return ("IGNORE:protonation_only_inchi (same connectivity "
                        "and stereo by InChI layer comparison, differs "
                        "only in protonation)")
            return ("IGNORE:stereo_only_inchi (same connectivity by InChI "
                    "layer comparison, differs in stereo layers)")

    # Same formula (or formulas match after normalization) — classify further
    s_mol = Chem.MolFromSmiles(stored_smiles) if stored_smiles else None
    p_mol = Chem.MolFromSmiles(pub_smiles) if pub_smiles else None
    if not s_mol or not p_mol:
        return "REVIEW:cannot_parse_smiles"

    # Check ring-chain tautomerism (e.g. cyclic sugar ↔ open-chain)
    s_ri = s_mol.GetRingInfo()
    p_ri = p_mol.GetRingInfo()
    s_o_rings = sum(
        1 for ring in s_ri.AtomRings()
        if any(s_mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring))
    p_o_rings = sum(
        1 for ring in p_ri.AtomRings()
        if any(p_mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring))
    if s_o_rings != p_o_rings:
        return (f"IGNORE:ring_chain_tautomer "
                f"(O-rings {s_o_rings}->{p_o_rings}, stored form is "
                f"biologically correct)")

    # Check tautomerism via RDKit canonical tautomer
    try:
        s_canon = Chem.MolToSmiles(_taut_enum.Canonicalize(s_mol))
        p_canon = Chem.MolToSmiles(_taut_enum.Canonicalize(p_mol))
        if s_canon == p_canon:
            return ("IGNORE:tautomer (same canonical tautomer, stored form "
                    "curated for biological relevance)")
    except Exception:
        pass

    # Check structural similarity to distinguish near-tautomers from
    # genuinely different compounds
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    fp1 = gen.GetFingerprint(s_mol)
    fp2 = gen.GetFingerprint(p_mol)
    sim = DataStructs.TanimotoSimilarity(fp1, fp2)

    # High similarity (>= TANIMOTO_CUTOFF) with same formula: likely a
    # tautomer or resonance form that RDKit's enumerator missed (common in
    # aromatic keto-enol systems like flavonoids)
    if sim >= TANIMOTO_CUTOFF:
        return (f"IGNORE:likely_tautomer (tanimoto={sim:.2f}, same formula, "
                f"high structural similarity — likely resonance/tautomer "
                f"not caught by RDKit enumerator)")

    # Low-to-medium similarity: genuinely different compound with same
    # formula (positional isomers, constitutional isomers)
    return (f"REVIEW:different_compound (tanimoto={sim:.2f}, same formula "
            f"but different connectivity — likely wrong xref mapping)")


def validate_correction(cpd_id, stored_struct, pubchem_struct, result_type):
    """Decide whether a PubChem correction should be accepted.

    Returns (accept: bool, reason: str).

    Rules:
    - MISMATCH: NEVER auto-corrected. Different InChIKey connectivity block
      means different molecular graph — cannot determine which is correct
      without domain expertise. Classified for manual review.
    - STEREO_DIFF: accepted only if PubChem has >= defined stereocenters
      (adds missing stereochemistry without losing existing assignments)
    - PROTONATION_DIFF: never auto-corrected (pH-dependent)
    - All: reject if PubChem SMILES is inconsistent with PubChem InChIKey
    """
    if result_type == "MISMATCH":
        detail = _classify_mismatch(stored_struct, pubchem_struct)
        return False, f"MISMATCH: {detail} — requires manual review"

    if result_type != "STEREO_DIFF":
        return False, f"Not auto-correctable ({result_type})"

    # Cross-check: verify PubChem SMILES and InChIKey describe the same molecule
    pub_smiles = pubchem_struct.get("smiles", "")
    pub_inchikey = pubchem_struct.get("inchikey", "")
    consistent, detail = _verify_smiles_inchikey_consistency(
        pub_smiles, pub_inchikey)
    if not consistent:
        return False, (f"STEREO_DIFF rejected: PubChem SMILES/InChIKey "
                       f"inconsistent ({detail})")

    # STEREO_DIFF: accept only if PubChem has >= defined stereocenters
    stored_smiles = stored_struct.get("smiles", "")
    old_spec, old_pot = count_defined_stereo(stored_smiles)
    new_spec, new_pot = count_defined_stereo(pub_smiles)

    if old_spec is None or new_spec is None:
        return False, ("STEREO_DIFF rejected: cannot parse SMILES for "
                       "stereo comparison")
    if new_spec < old_spec:
        return False, (f"STEREO_DIFF rejected: PubChem has fewer defined "
                       f"stereocenters ({new_spec} vs {old_spec})")

    # Verify shared stereocenters have the same configuration
    compatible, compat_detail = _check_stereo_compatibility(
        stored_struct, pubchem_struct)
    if not compatible:
        return False, (f"STEREO_DIFF rejected: {compat_detail}")

    return True, (f"STEREO_DIFF accepted: PubChem stereo "
                  f"{new_spec}>={old_spec} (of {new_pot} potential), "
                  f"{compat_detail}")


# ---------------------------------------------------------------------------
# Worker functions (module-level for multiprocessing)
# ---------------------------------------------------------------------------

def classify_mismatch_worker(args):
    """Worker function for parallel mismatch classification."""
    cpd_id, stored_struct, pub_struct, out_row = args
    detail = _classify_mismatch(stored_struct, pub_struct)
    return cpd_id, out_row, detail


def validate_correction_worker(args):
    """Worker function for parallel correction validation."""
    cpd_id, stored, corr, result_type = args
    accept, reason = validate_correction(cpd_id, stored, corr, result_type)
    return cpd_id, accept, reason
