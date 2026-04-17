#!/usr/bin/env python3
"""PubChem validation and correction pipeline for ModelSEED compound structures.

Validates compound structures against PubChem, classifies differences, applies
safe auto-corrections, and flags ambiguous cases for manual review.

Pipeline phases:
  Phase 0 — Self-consistency: validate/fix SMILES, InChI, InChIKey, formula, charge
  Phase 1 — External ID lookup: ChEBI/KEGG xref → PubChem CID (with conflict resolution)
  Phase 2 — Name lookup: batch name queries with disambiguation
  Phase 3 — Recovery: direct InChIKey/InChI lookup for unresolved compounds
  Phase 4 — Classification: classify results, sub-classify MISMATCHes,
             reclassify correctable sub-types, generate reports
  Phase 5 — Corrections: validate and apply STEREO_DIFF + PROTONATION_DIFF + pKa
  Phase 6 — Output: write both All and Unique format files, canonicalize SMILES

Usage:
  python Validate_PubChem_Structures.py                  # full run (resumable)
  python Validate_PubChem_Structures.py --limit 100      # test with first 100
  python Validate_PubChem_Structures.py --apply          # run + apply corrections
  python Validate_PubChem_Structures.py --test           # run unit tests
  python Validate_PubChem_Structures.py --clear-cache    # fresh re-validation
"""

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Imports, Constants, Path Resolution, CLI Args
# ═══════════════════════════════════════════════════════════════════════════════

import argparse
import csv
import glob
import logging
import os
import re
import sqlite3
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import Pool

import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import FindPotentialStereo, StereoSpecified
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey, MolFromInchi
from rdkit.Chem.MolStandardize import rdMolStandardize

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python'))
from BiochemPy import Compounds, InChIs

DB_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))
STRUCTURES_DIR = os.path.join(DB_ROOT, 'Structures')
ALIASES_DIR = os.path.join(DB_ROOT, 'Aliases')

ALL_STRUCTURES_FILE = os.path.join(STRUCTURES_DIR, 'All_ModelSEED_Structures.txt')
UNIQUE_STRUCTURES_FILE = os.path.join(STRUCTURES_DIR, 'Unique_ModelSEED_Structures.txt')
NAMES_FILE = os.path.join(ALIASES_DIR, 'Unique_ModelSEED_Compound_Names.txt')
ALIASES_FILE = os.path.join(ALIASES_DIR, 'Unique_ModelSEED_Compound_Aliases.txt')

CACHE_DB = os.path.join(SCRIPT_DIR, 'pubchem_cache.sqlite')
LOG_FILE = os.path.join(SCRIPT_DIR, 'pubchem_validate.log')
REPORT_FILE = os.path.join(SCRIPT_DIR, 'pubchem_validation_report.tsv')
MISMATCH_FILE = os.path.join(SCRIPT_DIR, 'pubchem_mismatches.tsv')
CORRECTIONS_LOG = os.path.join(SCRIPT_DIR, 'pubchem_corrections_log.tsv')
PROTONATION_FILE = os.path.join(SCRIPT_DIR, 'pubchem_protonation_diffs.tsv')
STEREO_FILE = os.path.join(SCRIPT_DIR, 'pubchem_stereo_diffs.tsv')
REJECTED_FILE = os.path.join(SCRIPT_DIR, 'pubchem_rejected_corrections.tsv')
PROTONATION_CORRECTIONS_FILE = os.path.join(SCRIPT_DIR,
                                             'pubchem_protonation_corrections.tsv')
REVIEW_FILE = os.path.join(SCRIPT_DIR, 'pubchem_review_different_compounds.tsv')
CONSISTENCY_FILE = os.path.join(SCRIPT_DIR, 'consistency_report.tsv')
TAUTOMER_FILE = os.path.join(SCRIPT_DIR, 'pubchem_tautomer_diffs.tsv')
IMAGES_DIR = os.path.join(SCRIPT_DIR, 'struct_imgs')

logger = logging.getLogger("pubchem_validate")
logger.setLevel(logging.DEBUG)
_log_formatter = logging.Formatter(
    "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_RETRIES = 4
NAME_BATCH_SIZE = 100
CID_BATCH_SIZE = 100
TANIMOTO_CUTOFF = 0.8
MAX_NAMES = 3

_RESULT_PRIORITY = {"MATCH": 0, "PROTONATION_DIFF": 1, "STEREO_DIFF": 2,
                     "MISMATCH": 3, "NO_KEY": 4}

_taut_enum = rdMolStandardize.TautomerEnumerator()
_uncharger = rdMolStandardize.Uncharger(canonicalOrder=False)
_morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)


def _parse_formula(f):
    return Counter({m[0]: int(m[1]) if m[1] else 1
                    for m in re.findall(r'([A-Z][a-z]?)(\d*)', f) if m[0]})


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate and correct ModelSEED compound structures "
                    "via PubChem.")
    parser.add_argument("--limit", type=int, default=0,
                        help="Only process first N compounds (0 = all)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from cache (default behavior)")
    parser.add_argument("--workers", type=int, default=10,
                        help="Number of concurrent request threads (default 10)")
    parser.add_argument("--apply", action="store_true",
                        help="Fetch full structures and apply corrections")
    parser.add_argument("--images", action="store_true",
                        help="Generate structure comparison images")
    parser.add_argument("--max-images", type=int, default=20,
                        help="Max images per result type (default 20)")
    parser.add_argument("--clear-cache", action="store_true",
                        help="Clear the SQLite cache for a fresh run")
    parser.add_argument("--skip-ph7", action="store_true",
                        help="Skip pH 7 protonation normalization")
    parser.add_argument("--review-report", action="store_true",
                        help="Generate detailed review for different_compound "
                             "mismatches")
    parser.add_argument("--test", action="store_true",
                        help="Run built-in unit tests and exit")
    return parser.parse_args()


def _write_tsv(path, header, rows, sort_key=None):
    if sort_key and rows:
        rows = sorted(rows, key=sort_key)
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: PubChem API
# ═══════════════════════════════════════════════════════════════════════════════

def _extract_smiles(props):
    return (props.get("IsomericSMILES")
            or props.get("CanonicalSMILES")
            or props.get("SMILES", ""))


class RateLimiter:
    def __init__(self, rate=5.0):
        self._interval = 1.0 / rate
        self._lock = threading.Lock()
        self._next_time = time.monotonic()

    def acquire(self):
        with self._lock:
            now = time.monotonic()
            if now < self._next_time:
                wait = self._next_time - now
                self._next_time += self._interval
            else:
                wait = 0
                self._next_time = now + self._interval
        if wait > 0:
            time.sleep(wait)


_rate_limiter = RateLimiter(rate=5.0)


def pubchem_request(method, url, **kwargs):
    for attempt in range(MAX_RETRIES):
        _rate_limiter.acquire()
        try:
            resp = (requests.get(url, timeout=30, **kwargs) if method == "GET"
                    else requests.post(url, timeout=30, **kwargs))
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                return None
            if resp.status_code in (429, 500, 502, 503, 504):
                wait = 2 ** (attempt + 1)
                logger.warning("PubChem %d, retrying in %ds...",
                               resp.status_code, wait)
                time.sleep(wait)
                continue
            return None
        except (requests.RequestException, ValueError):
            time.sleep(2 ** (attempt + 1))
    return None


def _pubchem_props(method, url, **kwargs):
    data = pubchem_request(method, url, **kwargs)
    if not data or "PropertyTable" not in data:
        return None
    props = data["PropertyTable"]["Properties"]
    if not props:
        return None
    p = props[0]
    return (str(p.get("CID", "")), p.get("InChIKey", ""), _extract_smiles(p))


def query_xref(xref_id):
    return _pubchem_props(
        "GET", f"{PUBCHEM_BASE}/compound/xref/RegistryID/"
               f"{xref_id}/property/InChIKey,IsomericSMILES/JSON")


def query_inchikey(inchikey):
    return _pubchem_props(
        "GET", f"{PUBCHEM_BASE}/compound/inchikey/{inchikey}"
               f"/property/InChIKey,IsomericSMILES/JSON")


def query_inchi(inchi_str):
    return _pubchem_props(
        "POST", f"{PUBCHEM_BASE}/compound/inchi"
                f"/property/InChIKey,IsomericSMILES/JSON",
        data={"inchi": inchi_str})


def _query_names_batch_recursive(name_list):
    if not name_list:
        return {}
    url = (f"{PUBCHEM_BASE}/compound/name"
           f"/property/InChIKey,IsomericSMILES/JSON")
    if len(name_list) == 1:
        name = name_list[0]
        data = pubchem_request("POST", url, data={"name": name})
        if data and "PropertyTable" in data:
            props = data["PropertyTable"]["Properties"]
            if not props:
                return {}
            p = props[0]
            return {name: (str(p.get("CID", "")),
                           p.get("InChIKey", ""), _extract_smiles(p))}
        return {}
    results = {}
    for name in name_list:
        r = _query_names_batch_recursive([name])
        results.update(r)
    return results


def query_cid_properties_batch(cid_list):
    if not cid_list:
        return {}
    url = (f"{PUBCHEM_BASE}/compound/cid"
           f"/property/InChIKey,IsomericSMILES,InChI/JSON")
    data = pubchem_request("POST", url,
                           data={"cid": ",".join(str(c) for c in cid_list)})
    results = {}
    if data and "PropertyTable" in data:
        for p in data["PropertyTable"]["Properties"]:
            cid = str(p.get("CID", ""))
            if cid:
                results[cid] = {"smiles": _extract_smiles(p),
                                "inchi": p.get("InChI", ""),
                                "inchikey": p.get("InChIKey", "")}
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: SQLite Cache
# ═══════════════════════════════════════════════════════════════════════════════

def init_cache(db_path):
    conn = sqlite3.connect(db_path, check_same_thread=False)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cache (
            cpd_id TEXT PRIMARY KEY, strategy TEXT, query TEXT, status TEXT,
            pubchem_cid TEXT, pubchem_inchikey TEXT, pubchem_smiles TEXT,
            timestamp REAL)""")
    conn.execute("""
        CREATE TABLE IF NOT EXISTS corrections (
            cpd_id TEXT PRIMARY KEY, pubchem_cid TEXT, pubchem_smiles TEXT,
            pubchem_inchi TEXT, pubchem_inchikey TEXT, timestamp REAL)""")
    conn.commit()
    return conn


def get_cached(conn):
    return {row[0] for row in conn.execute("SELECT cpd_id FROM cache")}


def save_batch_to_cache(conn, lock, rows):
    with lock:
        try:
            conn.executemany(
                "INSERT OR REPLACE INTO cache VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                rows)
            conn.commit()
        except sqlite3.Error:
            logger.warning("Batch cache write failed, falling back to "
                           "row-by-row for %d rows", len(rows))
            for row in rows:
                try:
                    conn.execute(
                        "INSERT OR REPLACE INTO cache "
                        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)", row)
                except sqlite3.Error:
                    logger.warning("Skipping bad cache row: %s",
                                   row[0] if row else "?")
            conn.commit()


# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Data Loading
# ═══════════════════════════════════════════════════════════════════════════════

def load_structures(path):
    structs = {}
    with open(path, "r") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 8 or row[2] != "Charged":
                continue
            cpd_id, typ = row[0], row[1]
            formula, charge, structure = row[5], row[6], row[7]
            if cpd_id not in structs:
                structs[cpd_id] = {}
            s = structs[cpd_id]
            if formula != "null" and "formula" not in s:
                s["formula"] = formula
                s["charge"] = charge
            key = typ.strip().lower()
            if structure and structure != "null":
                field = {"smile": "smiles", "inchikey": "inchikey",
                         "inchi": "inchi"}.get(key)
                if field and field not in s:
                    s[field] = structure
    return structs


def load_names(path):
    names = defaultdict(list)
    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) >= 2 and row[1].strip():
                names[row[0]].append(row[1].strip())
    return dict(names)


def load_external_ids(path):
    ids = defaultdict(lambda: {"chebi": [], "kegg": []})
    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 3:
                continue
            if row[2] == "ChEBI":
                ids[row[0]]["chebi"].append(row[1])
            elif row[2] == "KEGG":
                ids[row[0]]["kegg"].append(row[1])
    return dict(ids)


def load_pka_from_db():
    pka_data = {}
    for filepath in sorted(glob.glob(os.path.join(DB_ROOT, "compound_*.tsv"))):
        with open(filepath) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                cpd_id = row['id']
                charge = row.get('charge', '0')
                pka_data[cpd_id] = {
                    'pka': _parse_pka_string(row.get('pka', '') or ''),
                    'pkb': _parse_pka_string(row.get('pkb', '') or ''),
                    'db_charge': int(charge) if charge else 0}
    return pka_data


def _parse_pka_string(s):
    if not s or s == 'null':
        return []
    entries = []
    for part in s.split(';'):
        fields = part.split(':')
        if len(fields) == 3:
            try:
                entries.append(
                    (int(fields[0]), int(fields[1]), float(fields[2])))
            except ValueError:
                continue
    return entries


# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: Structure Comparison & Classification
# ═══════════════════════════════════════════════════════════════════════════════

def pick_best_names(name_list, max_names=MAX_NAMES):
    cleaned = [re.sub(r"<[^>]+>", "", n).strip() for n in name_list]
    cleaned = [n for n in cleaned if n]
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
    seen, result = set(), []
    for _, n in scored:
        key = n.lower()
        if key not in seen:
            seen.add(key)
            result.append(n)
            if len(result) >= max_names:
                break
    return result


def compare_inchikeys(stored, pubchem):
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
    if not stored_inchi or not pubchem_inchi:
        return None
    try:
        s_formula, s_layers = InChIs.parse(stored_inchi)
        p_formula, p_layers = InChIs.parse(pubchem_inchi)
    except Exception:
        return None
    s_no_stereo = InChIs.build(s_formula, s_layers,
                                remove=('b', 't', 'm', 's'))
    p_no_stereo = InChIs.build(p_formula, p_layers,
                                remove=('b', 't', 'm', 's'))
    s_no_prot = InChIs.build(s_formula, s_layers, remove=('p', 'q'))
    p_no_prot = InChIs.build(p_formula, p_layers, remove=('p', 'q'))
    s_conn = InChIs.build(s_formula, s_layers,
                           remove=('b', 't', 'm', 's', 'p', 'q'))
    p_conn = InChIs.build(p_formula, p_layers,
                           remove=('b', 't', 'm', 's', 'p', 'q'))
    return {
        "connectivity_match": s_conn == p_conn and s_conn != "",
        "stereo_match": s_no_prot == p_no_prot and s_no_prot != "",
        "protonation_match": s_no_stereo == p_no_stereo and s_no_stereo != "",
    }


def count_defined_stereo(smiles):
    if not smiles:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    stereo_info = FindPotentialStereo(mol)
    specified = sum(1 for s in stereo_info
                    if s.specified == StereoSpecified.Specified)
    return specified, len(stereo_info)


def _check_stereo_compatibility(stored_struct, pubchem_struct):
    stored_inchi = stored_struct.get("inchi", "")
    pub_inchi = pubchem_struct.get("inchi", "")
    if not stored_inchi or not pub_inchi:
        return True, "skipped (no InChI data)"
    try:
        _, s_layers = InChIs.parse(stored_inchi)
        _, p_layers = InChIs.parse(pub_inchi)
    except Exception:
        return True, "skipped (InChI parse error)"
    s_t, p_t = s_layers.get('t', ''), p_layers.get('t', '')
    s_b, p_b = s_layers.get('b', ''), p_layers.get('b', '')
    if not s_t and not s_b:
        return True, "compatible (stored has no InChI stereo layers)"
    inversions, checked = 0, 0
    if s_t and p_t:
        s_centers = {m.group(1): m.group(2)
                     for m in re.finditer(r'(\d+)([+-])', s_t)}
        p_centers = {m.group(1): m.group(2)
                     for m in re.finditer(r'(\d+)([+-])', p_t)}
        for atom_num, s_config in s_centers.items():
            if atom_num in p_centers:
                checked += 1
                if s_config != p_centers[atom_num]:
                    inversions += 1
    if s_b and p_b:
        s_bonds = {m.group(1): m.group(2)
                   for m in re.finditer(r'(\d+-\d+)([+-])', s_b)}
        p_bonds = {m.group(1): m.group(2)
                   for m in re.finditer(r'(\d+-\d+)([+-])', p_b)}
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
    if not inchi_str:
        return None, None
    try:
        formula, layers = InChIs.parse(inchi_str)
        chg = InChIs.charge(layers['q'], layers['p'])
        if '.' in formula:
            formula, _ = Compounds.mergeFormula(formula)
        if layers['p']:
            formula, _ = InChIs.adjust_protons(formula, layers['p'])
        return formula, chg
    except Exception:
        return None, None


def _classify_mismatch(stored_struct, pubchem_struct):
    """Returns (detail_string, reclassify_as) where reclassify_as is
    'STEREO_DIFF', 'PROTONATION_DIFF', or None."""
    stored_smiles = stored_struct.get("smiles", "")
    pub_smiles = pubchem_struct.get("smiles", "")
    stored_formula = stored_struct.get("formula", "")
    pub_inchi = pubchem_struct.get("inchi", "")

    pub_formula, _ = compute_formula_charge_from_inchi(pub_inchi)
    if not pub_formula:
        pub_formula = pubchem_struct.get("formula", "")
    sf = (stored_formula or "").strip()
    pf = (pub_formula or "").strip()

    if sf and pf and sf != pf:
        sf_merged, _ = Compounds.mergeFormula(sf)
        pf_merged, _ = Compounds.mergeFormula(pf)
        if (sf_merged != "null" and pf_merged != "null"
                and sf_merged != pf_merged):
            sf_atoms = Compounds.parseFormula(sf_merged)
            pf_atoms = Compounds.parseFormula(pf_merged)
            sf_no_h = {k: v for k, v in sf_atoms.items() if k != "H"}
            pf_no_h = {k: v for k, v in pf_atoms.items() if k != "H"}
            if sf_no_h == pf_no_h:
                h_diff = pf_atoms.get("H", 0) - sf_atoms.get("H", 0)
                return (f"IGNORE:protonation_formula_diff (formula: "
                        f"stored={sf}, pubchem={pf}, H_diff={h_diff:+d})",
                        "PROTONATION_DIFF")
            return (f"IGNORE:wrong_mapping (formula: stored={sf}, "
                    f"pubchem={pf})", None)

    stored_inchi = stored_struct.get("inchi", "")
    if stored_inchi and pub_inchi:
        layer_cmp = compare_inchi_layers(stored_inchi, pub_inchi)
        if layer_cmp and layer_cmp["connectivity_match"]:
            if layer_cmp["stereo_match"]:
                return ("IGNORE:protonation_only_inchi (same connectivity "
                        "and stereo, differs in protonation)",
                        "PROTONATION_DIFF")
            return ("IGNORE:stereo_only_inchi (same connectivity, "
                    "differs in stereo layers)", "STEREO_DIFF")

    s_mol = Chem.MolFromSmiles(stored_smiles) if stored_smiles else None
    p_mol = Chem.MolFromSmiles(pub_smiles) if pub_smiles else None
    if not s_mol or not p_mol:
        return "REVIEW:cannot_parse_smiles", None

    s_o_rings = sum(
        1 for ring in s_mol.GetRingInfo().AtomRings()
        if any(s_mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring))
    p_o_rings = sum(
        1 for ring in p_mol.GetRingInfo().AtomRings()
        if any(p_mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in ring))
    if s_o_rings != p_o_rings:
        return (f"IGNORE:ring_chain_tautomer (O-rings "
                f"{s_o_rings}->{p_o_rings})"), None

    try:
        if (Chem.MolToSmiles(_taut_enum.Canonicalize(s_mol))
                == Chem.MolToSmiles(_taut_enum.Canonicalize(p_mol))):
            return ("IGNORE:tautomer (same canonical tautomer, stored "
                    "form curated for biological relevance)"), None
    except Exception:
        pass

    sim = DataStructs.TanimotoSimilarity(
        _morgan_gen.GetFingerprint(s_mol), _morgan_gen.GetFingerprint(p_mol))
    if sim >= TANIMOTO_CUTOFF:
        return (f"IGNORE:likely_tautomer (tanimoto={sim:.2f})"), None
    return (f"REVIEW:different_compound (tanimoto={sim:.2f}, likely "
            f"wrong xref mapping)"), None


def validate_correction(cpd_id, stored_struct, pubchem_struct, result_type):
    if result_type == "MISMATCH":
        detail, _ = _classify_mismatch(stored_struct, pubchem_struct)
        return False, f"MISMATCH: {detail} — requires manual review"
    if result_type != "STEREO_DIFF":
        return False, f"Not auto-correctable ({result_type})"
    pub_smiles = pubchem_struct.get("smiles", "")
    pub_inchikey = pubchem_struct.get("inchikey", "")
    # Verify PubChem SMILES/InChIKey consistency
    if pub_smiles and pub_inchikey:
        mol = Chem.MolFromSmiles(pub_smiles)
        if mol:
            derived = Chem.inchi.MolToInchiKey(mol)
            if derived and derived.split("-")[0] != pub_inchikey.split("-")[0]:
                return False, ("STEREO_DIFF rejected: PubChem SMILES/InChIKey "
                               "inconsistent")
        else:
            return False, "STEREO_DIFF rejected: cannot parse PubChem SMILES"
    stored_smiles = stored_struct.get("smiles", "")
    old_spec, old_pot = count_defined_stereo(stored_smiles)
    new_spec, new_pot = count_defined_stereo(pub_smiles)
    if old_spec is None or new_spec is None:
        return False, ("STEREO_DIFF rejected: cannot parse SMILES for "
                       "stereo comparison")
    if new_spec < old_spec:
        return False, (f"STEREO_DIFF rejected: PubChem has fewer defined "
                       f"stereocenters ({new_spec} vs {old_spec})")
    compatible, detail = _check_stereo_compatibility(
        stored_struct, pubchem_struct)
    if not compatible:
        return False, f"STEREO_DIFF rejected: {detail}"
    return True, (f"STEREO_DIFF accepted: PubChem stereo "
                  f"{new_spec}>={old_spec} (of {new_pot} potential), "
                  f"{detail}")


def classify_mismatch_worker(args):
    cpd_id, stored_struct, pub_struct, out_row = args
    detail, reclassify = _classify_mismatch(stored_struct, pub_struct)
    return cpd_id, out_row, detail, reclassify


def validate_correction_worker(args):
    cpd_id, stored, corr, result_type = args
    accept, reason = validate_correction(cpd_id, stored, corr, result_type)
    return cpd_id, accept, reason


# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Protonation (pH 7 normalization, pKa logic)
# ═══════════════════════════════════════════════════════════════════════════════

ACIDIC_SMARTS = [
    (Chem.MolFromSmarts('[SX4](=O)(=O)[OX2H1]'), 'sulfonic_acid'),
    (Chem.MolFromSmarts('[PX4](=O)[OX2H1]'), 'phosphoric_acid'),
    (Chem.MolFromSmarts('[CX3](=O)[OX2H1]'), 'carboxylic_acid'),
    (Chem.MolFromSmarts('[SX3](=O)[OX2H1]'), 'sulfinic_acid'),
    (Chem.MolFromSmarts('[#16X2H1]'), 'thiol'),
    (Chem.MolFromSmarts('[OX2H1][#15]'), 'phosphate_oh_generic'),
]

BASIC_SMARTS = [
    (Chem.MolFromSmarts('[NX3;H2;!$(NC=O);!$(N=*);!$(n)]'), 'primary_amine'),
    (Chem.MolFromSmarts('[NX3;H1;!$(NC=O);!$(N=*);!$(n)]'),
     'secondary_amine'),
    (Chem.MolFromSmarts('[NX3;H0;!$(NC=O);!$(N=*);!$(n)]'),
     'tertiary_amine'),
]

_STRONG_ACID_SMARTS = [
    Chem.MolFromSmarts('[PX4](=O)'),
    Chem.MolFromSmarts('[SX4](=O)(=O)'),
]


def predict_charge_from_pka(pka_vals, pkb_vals, ph=7.0):
    deprotonations = sum(1 for _, _, v in pka_vals if v < ph)
    protonations = sum(1 for _, _, v in pkb_vals if v > ph)
    return deprotonations, protonations, -deprotonations + protonations


def _should_skip_correction(mol, pka_info, stored_charge, db_charge):
    pka_vals, pkb_vals = pka_info['pka'], pka_info['pkb']
    if pka_vals or pkb_vals:
        _, _, pka_ionic = predict_charge_from_pka(pka_vals, pkb_vals)
        if abs(stored_charge - pka_ionic) < abs(db_charge - pka_ionic):
            return True, (f"pKa cross-validation: pKa_ionic={pka_ionic} "
                          f"closer to stored={stored_charge} than "
                          f"db_charge={db_charge}")
    else:
        for pattern in _STRONG_ACID_SMARTS:
            if mol.HasSubstructMatch(pattern):
                return True, ("no pKa data but molecule has strongly "
                              "acidic groups")
    return False, None


def adjust_smiles_to_target_charge(smiles, target_charge):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if Chem.GetFormalCharge(mol) == target_charge:
        return _mol_to_result(mol)
    try:
        mol = _uncharger.uncharge(mol)
    except Exception:
        return None
    charge_needed = target_charge - Chem.GetFormalCharge(mol)
    if charge_needed < 0:
        mol = _deprotonate(mol, abs(charge_needed))
    elif charge_needed > 0:
        mol = _protonate(mol, charge_needed)
    if mol is None or Chem.GetFormalCharge(mol) != target_charge:
        return None
    return _mol_to_result(mol)


def _deprotonate(mol, count):
    emol = Chem.RWMol(mol)
    removed = 0
    for pattern, _name in ACIDIC_SMARTS:
        if removed >= count:
            break
        for match in emol.GetSubstructMatches(pattern):
            if removed >= count:
                break
            atom = emol.GetAtomWithIdx(match[-1])
            total_h = atom.GetTotalNumHs()
            if total_h > 0:
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(total_h - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                removed += 1
    if removed < count:
        return None
    try:
        Chem.SanitizeMol(emol)
        return emol.GetMol()
    except Exception:
        return None


def _protonate(mol, count):
    emol = Chem.RWMol(mol)
    added = 0
    for pattern, _name in BASIC_SMARTS:
        if added >= count:
            break
        for match in emol.GetSubstructMatches(pattern):
            if added >= count:
                break
            atom = emol.GetAtomWithIdx(match[0])
            if atom.GetFormalCharge() == 0:
                total_h = atom.GetTotalNumHs()
                atom.SetNoImplicit(True)
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(total_h + 1)
                added += 1
    if added < count:
        return None
    try:
        Chem.SanitizeMol(emol)
        return emol.GetMol()
    except Exception:
        return None


def _mol_to_result(mol):
    smiles = Chem.MolToSmiles(mol)
    inchi = Chem.MolToInchi(mol)
    if not inchi:
        return None
    inchikey = Chem.InchiToInchiKey(inchi)
    if not inchikey:
        return None
    return {'smiles': smiles, 'inchi': inchi, 'inchikey': inchikey}


# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Pipeline Phases 0–5
# ═══════════════════════════════════════════════════════════════════════════════

def run_phase0_consistency(structures, report_file=None):
    if report_file is None:
        report_file = CONSISTENCY_FILE
    logger.info("Phase 0: Self-consistency validation")
    stats = {'total': 0, 'smiles_parse_fail': 0, 'inchi_computed': 0,
             'inchikey_computed': 0, 'inchikey_mismatch': 0,
             'connectivity_mismatch': 0, 'formula_fixed': 0,
             'charge_fixed': 0, 'skipped_rgroup': 0}
    report_rows = []
    consistency_fixes = {}

    for cpd_id in sorted(structures.keys()):
        s = structures[cpd_id]
        stats['total'] += 1
        smiles = s.get('smiles', '')
        inchi = s.get('inchi', '')
        inchikey = s.get('inchikey', '')
        formula = s.get('formula', '')
        charge = s.get('charge', '')

        if not smiles:
            continue
        if '*' in smiles:
            stats['skipped_rgroup'] += 1
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            stats['smiles_parse_fail'] += 1
            report_rows.append({'cpd_id': cpd_id, 'issue': 'smiles_parse_fail',
                                'field': 'SMILES', 'stored': smiles[:80],
                                'computed': '', 'action': 'flagged'})
            continue

        cpd_fixes = {}

        if not inchi:
            try:
                computed_inchi = MolToInchi(mol)
                if computed_inchi:
                    s['inchi'] = inchi = computed_inchi
                    stats['inchi_computed'] += 1
                    cpd_fixes['inchi'] = ('', computed_inchi)
                    report_rows.append({
                        'cpd_id': cpd_id, 'issue': 'inchi_computed',
                        'field': 'InChI', 'stored': '',
                        'computed': computed_inchi[:80], 'action': 'fixed'})
            except Exception:
                pass

        if inchi and not inchikey:
            try:
                computed_ik = InchiToInchiKey(inchi)
                if computed_ik:
                    s['inchikey'] = inchikey = computed_ik
                    stats['inchikey_computed'] += 1
                    cpd_fixes['inchikey'] = ('', computed_ik)
                    report_rows.append({
                        'cpd_id': cpd_id, 'issue': 'inchikey_computed',
                        'field': 'InChIKey', 'stored': '',
                        'computed': computed_ik, 'action': 'fixed'})
            except Exception:
                pass

        if not inchi or not inchikey:
            if cpd_fixes:
                consistency_fixes[cpd_id] = cpd_fixes
            continue

        try:
            computed_ik = InchiToInchiKey(inchi)
            if computed_ik and computed_ik != inchikey:
                stats['inchikey_mismatch'] += 1
                report_rows.append({
                    'cpd_id': cpd_id, 'issue': 'inchikey_mismatch',
                    'field': 'InChIKey', 'stored': inchikey,
                    'computed': computed_ik, 'action': 'flagged'})
        except Exception:
            pass

        try:
            smi_inchi = MolToInchi(mol)
            if smi_inchi:
                smi_ik = InchiToInchiKey(smi_inchi)
                stored_conn = inchikey.split('-')[0]
                smi_conn = smi_ik.split('-')[0] if smi_ik else ''
                if smi_conn and stored_conn != smi_conn:
                    stats['connectivity_mismatch'] += 1
                    report_rows.append({
                        'cpd_id': cpd_id,
                        'issue': 'smiles_inchi_connectivity_mismatch',
                        'field': 'SMILES/InChI',
                        'stored': f'InChI_conn={stored_conn}',
                        'computed': f'SMILES_conn={smi_conn}',
                        'action': 'flagged'})
        except Exception:
            pass

        comp_formula, comp_charge = compute_formula_charge_from_inchi(inchi)
        if comp_formula:
            if comp_formula != formula:
                if _parse_formula(formula) != _parse_formula(comp_formula):
                    old_formula = formula
                    s['formula'] = comp_formula
                    stats['formula_fixed'] += 1
                    cpd_fixes['formula'] = (old_formula, comp_formula)
                    report_rows.append({
                        'cpd_id': cpd_id, 'issue': 'formula_mismatch',
                        'field': 'Formula', 'stored': old_formula,
                        'computed': comp_formula, 'action': 'fixed'})
            comp_charge_str = (str(comp_charge)
                               if comp_charge is not None else '')
            if comp_charge_str and comp_charge_str != str(charge):
                old_charge = charge
                s['charge'] = comp_charge_str
                stats['charge_fixed'] += 1
                cpd_fixes['charge'] = (str(old_charge), comp_charge_str)
                report_rows.append({
                    'cpd_id': cpd_id, 'issue': 'charge_mismatch',
                    'field': 'Charge', 'stored': str(old_charge),
                    'computed': comp_charge_str, 'action': 'fixed'})

        if cpd_fixes:
            consistency_fixes[cpd_id] = cpd_fixes

    if report_rows:
        fieldnames = ['cpd_id', 'issue', 'field', 'stored', 'computed',
                      'action']
        with open(report_file, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in sorted(report_rows, key=lambda r: r['cpd_id']):
                writer.writerow(row)
        logger.info("  Consistency report: %s", report_file)

    logger.info("  Total compounds: %d", stats['total'])
    logger.info("  Skipped (R-groups): %d", stats['skipped_rgroup'])
    for key in ('smiles_parse_fail', 'inchi_computed', 'inchikey_computed',
                'inchikey_mismatch', 'connectivity_mismatch',
                'formula_fixed', 'charge_fixed'):
        if stats[key]:
            logger.info("  %s: %d", key.replace('_', ' '), stats[key])
    logger.info("  Consistency fixes to persist: %d", len(consistency_fixes))
    return consistency_fixes


def run_phase1_pubchem_lookup(to_process, external_ids, structures, conn,
                              db_lock, workers, candidates):
    chebi_to_cpds = defaultdict(list)
    kegg_to_cpds = defaultdict(list)
    cpd_has_chebi, cpd_has_kegg = set(), set()
    needs_name_lookup = []

    for cpd_id in to_process:
        ext = external_ids.get(cpd_id, {"chebi": [], "kegg": []})
        has_chebi = has_kegg = False
        for cid_val in ext["chebi"]:
            chebi_to_cpds[f"CHEBI:{cid_val}"].append(cpd_id)
            has_chebi = True
        for kid in ext["kegg"]:
            kegg_to_cpds[kid].append(cpd_id)
            has_kegg = True
        if has_chebi:
            cpd_has_chebi.add(cpd_id)
        if has_kegg:
            cpd_has_kegg.add(cpd_id)
        if not has_chebi and not has_kegg:
            needs_name_lookup.append(cpd_id)

    unique_chebi = list(chebi_to_cpds.keys())
    unique_kegg = list(kegg_to_cpds.keys())
    total_xref_queries = len(unique_chebi) + len(unique_kegg)

    logger.info("Phase 1: External ID lookup (ChEBI + KEGG xrefs)")
    logger.info("  Unique ChEBI: %d, KEGG: %d, total queries: %d",
                len(unique_chebi), len(unique_kegg), total_xref_queries)

    xref_resolved_cpds = set()
    xref_conflicts = set()
    xref_conflict_pairs = defaultdict(list)
    pending_cache = []
    pending_index = {}
    CACHE_FLUSH_SIZE = 100

    def flush_cache():
        if pending_cache:
            save_batch_to_cache(conn, db_lock, list(pending_cache))
            pending_cache.clear()
            pending_index.clear()

    def do_xref_query(xref_id, strategy):
        return (xref_id, strategy, query_xref(xref_id))

    pbar = tqdm(total=total_xref_queries,
                desc="Phase 1: xref queries") if tqdm else None

    xref_done = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = []
        for cid in unique_chebi:
            futures.append(executor.submit(do_xref_query, cid, "chebi_xref"))
        for kid in unique_kegg:
            futures.append(executor.submit(do_xref_query, kid, "kegg_xref"))

        for future in as_completed(futures):
            try:
                xref_id, strategy, result = future.result()
            except Exception:
                logger.warning("xref query failed", exc_info=True)
                xref_done += 1
                if pbar:
                    pbar.update(1)
                continue
            xref_done += 1
            if pbar:
                pbar.update(1)
            if xref_done % 500 == 0:
                logger.info("  Phase 1 progress: %d/%d",
                            xref_done, total_xref_queries)
            if result is None:
                continue
            cid, inchikey, smiles = result
            cpd_ids = (chebi_to_cpds[xref_id] if strategy == "chebi_xref"
                       else kegg_to_cpds[xref_id])
            for cpd_id in cpd_ids:
                with db_lock:
                    existing = conn.execute(
                        "SELECT strategy, pubchem_cid FROM cache "
                        "WHERE cpd_id=?", (cpd_id,)).fetchone()
                    existing_full = None
                    if not existing and cpd_id in pending_index:
                        pi = pending_index[cpd_id]
                        existing = (pi[0], pi[1])
                        existing_full = pi
                if existing:
                    old_strategy, old_cid = existing
                    if (old_strategy != strategy and old_cid
                            and old_cid != cid):
                        xref_conflicts.add(cpd_id)
                        xref_resolved_cpds.discard(cpd_id)
                        if cpd_id not in xref_conflict_pairs:
                            if existing_full:
                                xref_conflict_pairs[cpd_id].append(
                                    (existing_full[1], existing_full[2],
                                     existing_full[3], existing_full[0],
                                     existing_full[4]))
                            else:
                                with db_lock:
                                    r = conn.execute(
                                        "SELECT pubchem_cid, "
                                        "pubchem_inchikey, pubchem_smiles, "
                                        "strategy, query FROM cache "
                                        "WHERE cpd_id=?",
                                        (cpd_id,)).fetchone()
                                if r:
                                    xref_conflict_pairs[cpd_id].append(
                                        (r[0], r[1], r[2], r[3], r[4]))
                        xref_conflict_pairs[cpd_id].append(
                            (cid, inchikey, smiles, strategy, xref_id))
                        query_str = (f"{old_strategy}->CID{old_cid};"
                                     f"{strategy}={xref_id}->CID{cid}")
                        with db_lock:
                            pending_cache.append((
                                cpd_id, "xref_conflict", query_str,
                                "xref_conflict", None, None, None,
                                time.time()))
                            pending_index[cpd_id] = (
                                "xref_conflict", None, None, None, None)
                else:
                    with db_lock:
                        pending_cache.append((
                            cpd_id, strategy, xref_id, "found",
                            cid, inchikey, smiles, time.time()))
                        pending_index[cpd_id] = (
                            strategy, cid, inchikey, smiles, xref_id)
                    xref_resolved_cpds.add(cpd_id)
                if len(pending_cache) >= CACHE_FLUSH_SIZE:
                    flush_cache()

    flush_cache()
    if pbar:
        pbar.close()

    logger.info("  Resolved via xref: %d", len(xref_resolved_cpds))
    if xref_conflicts:
        logger.info("  xref conflicts: %d", len(xref_conflicts))

    # Resolve conflicts inline
    if xref_conflict_pairs:
        resolution_rows = []
        resolved_count = 0
        for cpd_id, conflict_candidates in xref_conflict_pairs.items():
            stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
            if not stored_ik:
                continue
            scored = []
            for c_cid, c_ik, c_smi, c_strat, c_xref in conflict_candidates:
                if not c_ik:
                    continue
                rt = compare_inchikeys(stored_ik, c_ik)
                pri = _RESULT_PRIORITY.get(rt, 99)
                src = 0 if "chebi" in c_strat else 1
                scored.append((pri, src, c_cid, c_ik, c_smi,
                               c_strat, c_xref, rt))
            if not scored:
                continue
            scored.sort()
            best = scored[0]
            if best[0] <= 2:
                resolution_rows.append((
                    cpd_id, "xref_resolved",
                    f"xref_resolved:{best[5]}={best[6]}->CID{best[2]}"
                    f"({best[7]})",
                    "found", best[2], best[3], best[4], time.time()))
                xref_resolved_cpds.add(cpd_id)
                xref_conflicts.discard(cpd_id)
                resolved_count += 1
        if resolution_rows:
            save_batch_to_cache(conn, db_lock, resolution_rows)
        logger.info("  Conflicts resolved (InChIKey scoring): %d",
                    resolved_count)

    _resolve_cached_xref_conflicts(
        conn, db_lock, candidates, structures, workers)

    all_xref_cpds = set()
    for cpd_ids in chebi_to_cpds.values():
        all_xref_cpds.update(cpd_ids)
    for cpd_ids in kegg_to_cpds.values():
        all_xref_cpds.update(cpd_ids)
    unresolved_xref = [c for c in to_process if c in all_xref_cpds
                       and c not in xref_resolved_cpds
                       and c not in xref_conflicts]
    needs_name_lookup.extend(unresolved_xref)
    logger.info("  Compounds needing name lookup: %d", len(needs_name_lookup))
    return needs_name_lookup


def _resolve_cached_xref_conflicts(conn, db_lock, candidates, structures,
                                    workers):
    cur = conn.execute(
        "SELECT cpd_id, query FROM cache WHERE status='xref_conflict' "
        "AND cpd_id IN ({})".format(",".join("?" * len(candidates))),
        candidates)
    conflicts = cur.fetchall()
    if not conflicts:
        return
    logger.info("  Resolving %d cached xref conflicts via CID properties",
                len(conflicts))

    cpd_cids = {}
    all_cids = set()
    for cpd_id, query_str in conflicts:
        pairs = []
        for cid_str in re.findall(r'CID(\d+)', query_str):
            all_cids.add(cid_str)
            idx = query_str.index(f'CID{cid_str}')
            source = ('chebi' if 'chebi' in query_str[:idx].split(';')[-1]
                      else 'kegg')
            pairs.append((cid_str, source))
        cpd_cids[cpd_id] = pairs

    cid_props = {}
    batches = [sorted(all_cids)[i:i + CID_BATCH_SIZE]
               for i in range(0, len(all_cids), CID_BATCH_SIZE)]

    pbar = tqdm(total=len(batches),
                desc="  CID property batches") if tqdm else None
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(query_cid_properties_batch, b): b
                   for b in batches}
        for future in as_completed(futures):
            try:
                cid_props.update(future.result())
            except Exception:
                logger.warning("CID batch query failed", exc_info=True)
            if pbar:
                pbar.update(1)
    if pbar:
        pbar.close()

    resolution_rows = []
    resolved_count = 0
    for cpd_id, query_str in conflicts:
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        if not stored_ik:
            continue
        scored = []
        for cid_str, source in cpd_cids[cpd_id]:
            props = cid_props.get(cid_str)
            if not props or not props.get("inchikey"):
                continue
            rt = compare_inchikeys(stored_ik, props["inchikey"])
            pri = _RESULT_PRIORITY.get(rt, 99)
            scored.append((pri, 0 if source == 'chebi' else 1,
                           cid_str, props, rt))
        if not scored:
            continue
        scored.sort()
        best_pri, _, best_cid, best_props, best_rt = scored[0]
        if best_pri <= 2:
            resolution_rows.append((
                cpd_id, "xref_resolved",
                f"xref_resolved:CID{best_cid}({best_rt})", "found",
                best_cid, best_props["inchikey"],
                best_props.get("smiles", ""), time.time()))
            resolved_count += 1
    if resolution_rows:
        save_batch_to_cache(conn, db_lock, resolution_rows)
    logger.info("  Resolved via CID properties: %d", resolved_count)


def run_phase2_name_lookup(needs_name_lookup, names, structures, conn,
                           db_lock, workers):
    name_to_cpds = defaultdict(list)
    cpd_to_names = {}
    no_name = []

    for cpd_id in needs_name_lookup:
        best = pick_best_names(names.get(cpd_id, []))
        if best:
            cpd_to_names[cpd_id] = best
            for name in best:
                name_to_cpds[name].append(cpd_id)
        else:
            no_name.append(cpd_id)

    name_list = list(name_to_cpds.keys())
    logger.info("Phase 2: Name lookup")
    logger.info("  Unique names: %d, no usable name: %d",
                len(name_list), len(no_name))

    name_results = {}
    batches = [name_list[i:i + NAME_BATCH_SIZE]
               for i in range(0, len(name_list), NAME_BATCH_SIZE)]
    pbar = tqdm(total=len(batches),
                desc="Phase 2: name batches") if tqdm else None
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_query_names_batch_recursive, b): b
                   for b in batches}
        for future in as_completed(futures):
            try:
                name_results.update(future.result())
            except Exception:
                logger.warning("Name batch query failed", exc_info=True)
            if pbar:
                pbar.update(1)
    if pbar:
        pbar.close()

    found_names = not_found_names = ambiguous_names = 0
    name_cache_rows = []

    for cpd_id, tried_names in cpd_to_names.items():
        found_cids = {}
        for name in tried_names:
            if name in name_results:
                cid, ik, smi = name_results[name]
                found_cids[cid] = (name, ik, smi)
        if not found_cids:
            name_cache_rows.append((cpd_id, "name_lookup", tried_names[0],
                                    "not_found", None, None, None,
                                    time.time()))
            not_found_names += 1
        elif len(found_cids) == 1:
            cid = list(found_cids.keys())[0]
            name, ik, smi = found_cids[cid]
            name_cache_rows.append((cpd_id, "name_lookup", name, "found",
                                    cid, ik, smi, time.time()))
            found_names += 1
        else:
            stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
            resolved = (stored_ik and _disambiguate_name_results(
                cpd_id, stored_ik, found_cids, structures, name_cache_rows))
            if resolved:
                found_names += 1
            else:
                query_str = ";".join(f"{n}->CID{c}"
                                     for c, (n, _, _) in found_cids.items())
                name_cache_rows.append(
                    (cpd_id, "name_lookup_ambiguous", query_str,
                     "ambiguous", None, None, None, time.time()))
                ambiguous_names += 1
        if len(name_cache_rows) >= 100:
            save_batch_to_cache(conn, db_lock, name_cache_rows)
            name_cache_rows = []

    for cpd_id in no_name:
        name_cache_rows.append((cpd_id, "name_lookup", "", "not_found",
                                None, None, None, time.time()))
        not_found_names += 1

    if name_cache_rows:
        save_batch_to_cache(conn, db_lock, name_cache_rows)
    logger.info("  Found: %d, not found: %d, ambiguous: %d",
                found_names, not_found_names, ambiguous_names)


def _disambiguate_name_results(cpd_id, stored_ik, found_cids, structures,
                               cache_rows):
    s_parts = stored_ik.split("-")
    for cid, (name, ik, smi) in found_cids.items():
        if ik and ik == stored_ik:
            cache_rows.append((cpd_id, "name_lookup", name, "found",
                               cid, ik, smi, time.time()))
            return True
    for cid, (name, ik, smi) in found_cids.items():
        if not ik:
            continue
        p_parts = ik.split("-")
        if (len(s_parts) >= 2 and len(p_parts) >= 2
                and s_parts[0] == p_parts[0] and s_parts[1] == p_parts[1]):
            cache_rows.append((cpd_id, "name_lookup", name, "found",
                               cid, ik, smi, time.time()))
            return True
    block1_matches = [(cid, name, ik, smi)
                      for cid, (name, ik, smi) in found_cids.items()
                      if ik and ik.split("-")[0] == s_parts[0]]
    if len(block1_matches) == 1:
        cid, name, ik, smi = block1_matches[0]
        cache_rows.append((cpd_id, "name_lookup", name, "found",
                           cid, ik, smi, time.time()))
        return True
    if len(block1_matches) > 1:
        stored_smi = structures.get(cpd_id, {}).get("smiles", "")
        s_mol = Chem.MolFromSmiles(stored_smi) if stored_smi else None
        if s_mol:
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
            s_fp = gen.GetFingerprint(s_mol)
            best_sim, best_match = -1, None
            for cid, name, ik, smi in block1_matches:
                p_mol = Chem.MolFromSmiles(smi) if smi else None
                if p_mol:
                    sim = DataStructs.TanimotoSimilarity(
                        s_fp, gen.GetFingerprint(p_mol))
                    if sim > best_sim:
                        best_sim, best_match = sim, (cid, name, ik, smi)
            if best_match:
                cid, name, ik, smi = best_match
                cache_rows.append((cpd_id, "name_lookup", name, "found",
                                   cid, ik, smi, time.time()))
                return True
    return False


def run_phase3_recovery(candidates, structures, conn, db_lock, workers):
    cur = conn.execute(
        "SELECT cpd_id, strategy, status, pubchem_cid, pubchem_inchikey "
        "FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))), candidates)
    ik_candidates = []
    for row in cur.fetchall():
        cpd_id, strategy, status, pub_cid, pub_ik = row
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        if not stored_ik:
            continue
        if status == "not_found":
            ik_candidates.append((cpd_id, stored_ik, "NOT_FOUND"))
        elif status == "found" and pub_ik:
            rt = compare_inchikeys(stored_ik, pub_ik)
            if rt in ("MISMATCH", "STEREO_DIFF"):
                ik_candidates.append((cpd_id, stored_ik, rt))

    ik_to_cpds = defaultdict(list)
    for cpd_id, stored_ik, old_rt in ik_candidates:
        ik_to_cpds[stored_ik].append((cpd_id, old_rt))
    unique_iks = list(ik_to_cpds.keys())

    logger.info("Phase 3: Recovery (InChIKey + InChI direct lookup)")
    logger.info("  InChIKey candidates: %d (%d unique)",
                len(ik_candidates), len(unique_iks))

    if unique_iks:
        ik_results = {}
        pbar = tqdm(total=len(unique_iks),
                    desc="Phase 3a: InChIKey lookups") if tqdm else None
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(
                lambda ik=ik: (ik, query_inchikey(ik))): ik
                for ik in unique_iks}
            for future in as_completed(futures):
                try:
                    ik, result = future.result()
                except Exception:
                    logger.warning("InChIKey lookup failed", exc_info=True)
                    if pbar:
                        pbar.update(1)
                    continue
                ik_results[ik] = result
                if pbar:
                    pbar.update(1)
        if pbar:
            pbar.close()

        ik_cache_rows = []
        recovered = defaultdict(int)
        for stored_ik, cpd_list in ik_to_cpds.items():
            result = ik_results.get(stored_ik)
            if result is None:
                continue
            cid, pub_ik, smiles = result
            if not pub_ik:
                continue
            for cpd_id, old_rt in cpd_list:
                new_rt = compare_inchikeys(stored_ik, pub_ik)
                improved = (old_rt == "NOT_FOUND" or
                            (old_rt == "MISMATCH" and new_rt in
                             ("MATCH", "PROTONATION_DIFF", "STEREO_DIFF")) or
                            (old_rt == "STEREO_DIFF" and new_rt in
                             ("MATCH", "PROTONATION_DIFF")))
                if improved:
                    ik_cache_rows.append(
                        (cpd_id, "inchikey_lookup", stored_ik, "found",
                         cid, pub_ik, smiles, time.time()))
                    recovered[new_rt] += 1
        if ik_cache_rows:
            save_batch_to_cache(conn, db_lock, ik_cache_rows)
        logger.info("  InChIKey recovered: %d", sum(recovered.values()))

    # Part B: InChI recovery for still-not-found
    cur2 = conn.execute(
        "SELECT cpd_id, status FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))), candidates)
    inchi_candidates = [(cpd_id, structures.get(cpd_id, {}).get("inchi", ""))
                        for cpd_id, status in cur2.fetchall()
                        if status == "not_found"
                        and structures.get(cpd_id, {}).get("inchi", "")]
    logger.info("  InChI recovery candidates: %d", len(inchi_candidates))

    if inchi_candidates:
        inchi_to_cpds = defaultdict(list)
        for cpd_id, stored_inchi in inchi_candidates:
            inchi_to_cpds[stored_inchi].append(cpd_id)
        unique_inchis = list(inchi_to_cpds.keys())

        pbar = tqdm(total=len(unique_inchis),
                    desc="Phase 3b: InChI lookups") if tqdm else None
        inchi_results = {}
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(
                lambda i=i: (i, query_inchi(i))): i
                for i in unique_inchis}
            for future in as_completed(futures):
                try:
                    inchi_str, result = future.result()
                except Exception:
                    logger.warning("InChI lookup failed", exc_info=True)
                    if pbar:
                        pbar.update(1)
                    continue
                inchi_results[inchi_str] = result
                if pbar:
                    pbar.update(1)
        if pbar:
            pbar.close()

        inchi_cache_rows = []
        for stored_inchi, cpd_list in inchi_to_cpds.items():
            result = inchi_results.get(stored_inchi)
            if result is None:
                continue
            cid, pub_ik, smiles = result
            if not pub_ik:
                continue
            for cpd_id in cpd_list:
                inchi_cache_rows.append(
                    (cpd_id, "inchi_lookup", stored_inchi[:50], "found",
                     cid, pub_ik, smiles, time.time()))
        if inchi_cache_rows:
            save_batch_to_cache(conn, db_lock, inchi_cache_rows)
        logger.info("  InChI recovered: %d", len(inchi_cache_rows))


def run_phase5_corrections(conn, db_lock, candidates, structures, pka_data,
                           names, workers, skip_ph7=False,
                           reclassifications=None):
    if reclassifications is None:
        reclassifications = {}
    cur = conn.execute(
        "SELECT cpd_id, strategy, query, status, pubchem_cid, "
        "pubchem_inchikey FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))), candidates)
    correctable = {}
    protonation_diff_cpds = set()
    for row in cur.fetchall():
        cpd_id, strategy, query, status, pub_cid, pub_ik = row
        if status != "found" or not pub_cid:
            continue
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        rt = reclassifications.get(cpd_id) or compare_inchikeys(
            stored_ik, pub_ik)
        if rt in ("MISMATCH", "STEREO_DIFF"):
            correctable[cpd_id] = (pub_cid, strategy, query, rt)
        elif rt == "PROTONATION_DIFF":
            protonation_diff_cpds.add(cpd_id)

    logger.info("Phase 5: Corrections")
    logger.info("  STEREO_DIFF/MISMATCH candidates: %d", len(correctable))
    logger.info("  PROTONATION_DIFF candidates for pKa: %d",
                len(protonation_diff_cpds))

    # Batch-fetch structures for correctable compounds
    if correctable:
        already_corrected = {row[0] for row in
                             conn.execute("SELECT cpd_id FROM corrections")}
        unique_cids = {}
        for cpd_id, (cid, _, _, _) in correctable.items():
            if cpd_id not in already_corrected:
                unique_cids.setdefault(cid, []).append(cpd_id)

        if unique_cids:
            logger.info("  Fetching %d unique CID structures (batched)",
                        len(unique_cids))
            cid_list = sorted(unique_cids.keys())
            batches = [cid_list[i:i + CID_BATCH_SIZE]
                       for i in range(0, len(cid_list), CID_BATCH_SIZE)]
            pbar = tqdm(total=len(batches),
                        desc="Phase 5: CID batches") if tqdm else None
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {executor.submit(
                    query_cid_properties_batch, b): b for b in batches}
                for future in as_completed(futures):
                    try:
                        batch_results = future.result()
                    except Exception:
                        logger.warning("CID property batch failed",
                                       exc_info=True)
                        if pbar:
                            pbar.update(1)
                        continue
                    if pbar:
                        pbar.update(1)
                    rows = []
                    for cid, props in batch_results.items():
                        for cpd_id in unique_cids.get(cid, []):
                            rows.append((cpd_id, cid, props["smiles"],
                                         props["inchi"], props["inchikey"],
                                         time.time()))
                    if rows:
                        with db_lock:
                            conn.executemany(
                                "INSERT OR REPLACE INTO corrections "
                                "VALUES (?,?,?,?,?,?)", rows)
                            conn.commit()
            if pbar:
                pbar.close()

    # Build corrections from DB
    raw_result = {}
    if correctable:
        cur3 = conn.execute(
            "SELECT cpd_id, pubchem_cid, pubchem_smiles, pubchem_inchi, "
            "pubchem_inchikey FROM corrections WHERE cpd_id IN ({})".format(
                ",".join("?" * len(correctable))), list(correctable.keys()))
        for row in cur3.fetchall():
            cpd_id = row[0]
            cid, strategy, query, rt = correctable[cpd_id]
            raw_result[cpd_id] = {
                "pubchem_cid": row[1], "smiles": row[2], "inchi": row[3],
                "inchikey": row[4], "result_type": rt,
                "strategy": strategy, "query": query}

    # Validate corrections
    work_items = [(cpd_id, structures.get(cpd_id, {}), corr,
                   corr["result_type"])
                  for cpd_id, corr in raw_result.items()]
    validated = {}
    rejected = []
    if work_items:
        num_workers = min(os.cpu_count() or 4, 32)
        with Pool(num_workers) as pool:
            for cpd_id, accept, reason in pool.imap_unordered(
                    validate_correction_worker, work_items, chunksize=32):
                corr = raw_result[cpd_id]
                if accept:
                    corr["validation_reason"] = reason
                    validated[cpd_id] = corr
                else:
                    rejected.append((cpd_id, corr["result_type"], reason))

    logger.info("  STEREO_DIFF validated: %d, rejected: %d",
                len(validated), len(rejected))

    if rejected:
        rows = []
        for cpd_id, rt, reason in rejected:
            corr = raw_result.get(cpd_id, {})
            stored = structures.get(cpd_id, {})
            cpd_names = names.get(cpd_id, []) if names else []
            rows.append([cpd_id, cpd_names[0] if cpd_names else "", rt,
                         reason, corr.get("pubchem_cid", ""),
                         corr.get("strategy", ""), corr.get("query", ""),
                         stored.get("smiles", ""), corr.get("smiles", ""),
                         stored.get("inchikey", ""),
                         corr.get("inchikey", "")])
        _write_tsv(REJECTED_FILE,
                   ["cpd_id", "compound_name", "result_type", "reason",
                    "pubchem_cid", "strategy", "query", "stored_smiles",
                    "pubchem_smiles", "stored_inchikey", "pubchem_inchikey"],
                   rows)

    # Normalize STEREO_DIFF corrections to pH 7
    if validated and not skip_ph7:
        normalized = 0
        for cpd_id, corr in list(validated.items()):
            if corr.get("result_type") != "STEREO_DIFF":
                continue
            pub_smiles = corr.get("smiles", "")
            pka_info = pka_data.get(cpd_id)
            if not pub_smiles or pka_info is None:
                continue
            target = pka_info.get('db_charge')
            if target is None:
                continue
            ph7 = adjust_smiles_to_target_charge(pub_smiles, target)
            if ph7 is None:
                continue
            if compare_inchikeys(corr.get("inchikey", ""),
                                 ph7["inchikey"]) == "PROTONATION_DIFF":
                corr.update(ph7)
                normalized += 1
        if normalized:
            logger.info("  STEREO_DIFF normalized to pH 7: %d", normalized)

    # PROTONATION_DIFF via pKa
    if not skip_ph7 and protonation_diff_cpds:
        prot_corrected = 0
        for cpd_id in protonation_diff_cpds:
            if cpd_id in validated:
                continue
            pka_info = pka_data.get(cpd_id)
            if pka_info is None:
                continue
            stored = structures.get(cpd_id, {})
            stored_smiles = stored.get("smiles", "")
            if not stored_smiles:
                continue
            mol = Chem.MolFromSmiles(stored_smiles)
            if mol is None:
                continue
            stored_charge = Chem.GetFormalCharge(mol)
            target_charge = pka_info.get('db_charge')
            if target_charge is None or stored_charge == target_charge:
                continue
            if any(6.0 <= v <= 8.0 for _, _, v in pka_info['pka']):
                continue
            skip, _ = _should_skip_correction(
                mol, pka_info, stored_charge, target_charge)
            if skip:
                continue
            result = adjust_smiles_to_target_charge(
                stored_smiles, target_charge)
            if result is None:
                continue
            if compare_inchikeys(stored.get("inchikey", ""),
                                 result['inchikey']) not in (
                    "PROTONATION_DIFF", "MATCH"):
                continue
            validated[cpd_id] = {
                'smiles': result['smiles'], 'inchi': result['inchi'],
                'inchikey': result['inchikey'],
                'pubchem_cid': 'pKa+PubChem',
                'result_type': 'PROTONATION_DIFF_CORRECTED',
                'strategy': 'pubchem_protonation_pka',
                'query': 'PubChem_PROTONATION_DIFF+ChemAxon_pKa',
                'validation_reason': (
                    f"PROTONATION_DIFF confirmed by pKa: "
                    f"stored={stored_charge}, target={target_charge}")}
            prot_corrected += 1
        logger.info("  PROTONATION_DIFF corrected via pKa: %d",
                    prot_corrected)

    # General pKa corrections
    if not skip_ph7:
        pka_corrections = _run_pka_validation(
            candidates, structures, pka_data, names)
        pka_added = sum(1 for c in pka_corrections if c not in validated)
        for cpd_id, corr in pka_corrections.items():
            if cpd_id not in validated:
                validated[cpd_id] = corr
        if pka_added:
            logger.info("  pKa-only corrections: %d", pka_added)

    return validated


def _run_pka_validation(candidates, structures, pka_data, names):
    work_items = [(cpd_id, structures.get(cpd_id, {}).get("smiles", ""))
                  for cpd_id in candidates
                  if structures.get(cpd_id, {}).get("smiles", "")]

    corrections = {}
    stats = {'already_correct': 0, 'pka_corrected': 0, 'no_pka_data': 0,
             'skipped_crossval': 0, 'adjustment_failed': 0,
             'parse_failed': 0, 'borderline_skipped': 0}
    correction_rows = []

    items_iter = (tqdm(work_items, desc="Phase 5: pKa validation")
                  if tqdm else work_items)

    for cpd_id, stored_smiles in items_iter:
        pka_info = pka_data.get(cpd_id)
        if pka_info is None:
            stats['no_pka_data'] += 1
            continue
        mol = Chem.MolFromSmiles(stored_smiles)
        if mol is None:
            stats['parse_failed'] += 1
            continue
        stored_charge = Chem.GetFormalCharge(mol)
        target_charge = pka_info.get('db_charge')
        if target_charge is None:
            stats['parse_failed'] += 1
            continue

        borderline = [(a, v) for _, a, v in pka_info['pka']
                      if 6.0 <= v <= 8.0]
        if borderline:
            if stored_charge != target_charge:
                stats['borderline_skipped'] += 1
                stored = structures.get(cpd_id, {})
                cpd_names = names.get(cpd_id, []) if names else []
                correction_rows.append({
                    'cpd_id': cpd_id,
                    'compound_name': cpd_names[0] if cpd_names else "",
                    'action': 'flagged_borderline',
                    'stored_inchikey': stored.get('inchikey', ''),
                    'pka_inchikey': '', 'stored_smiles': stored_smiles,
                    'pka_smiles': '',
                    'stored_formula': stored.get('formula', ''),
                    'pka_formula': '',
                    'stored_charge': str(stored_charge),
                    'pka_charge': str(target_charge),
                    'borderline_pka': "; ".join(
                        f"atom {a} pKa={v:.2f}" for a, v in borderline)})
            continue

        if stored_charge == target_charge:
            stats['already_correct'] += 1
            continue
        skip, _ = _should_skip_correction(
            mol, pka_info, stored_charge, target_charge)
        if skip:
            stats['skipped_crossval'] += 1
            continue
        result = adjust_smiles_to_target_charge(stored_smiles, target_charge)
        if result is None:
            stats['adjustment_failed'] += 1
            continue
        stored = structures.get(cpd_id, {})
        stored_ik = stored.get("inchikey", "")
        if compare_inchikeys(stored_ik, result['inchikey']) not in (
                "PROTONATION_DIFF", "MATCH"):
            stats['adjustment_failed'] += 1
            continue
        formula, _ = compute_formula_charge_from_inchi(result['inchi'])
        corrections[cpd_id] = {
            'smiles': result['smiles'], 'inchi': result['inchi'],
            'inchikey': result['inchikey'], 'pubchem_cid': 'pKa',
            'result_type': 'PKA_CORRECTION', 'strategy': 'chemaxon_pka',
            'query': 'ChemAxon_pKa',
            'validation_reason': (f"Protonation adjusted "
                                  f"(stored={stored_charge}, "
                                  f"target={target_charge})")}
        cpd_names = names.get(cpd_id, []) if names else []
        correction_rows.append({
            'cpd_id': cpd_id,
            'compound_name': cpd_names[0] if cpd_names else "",
            'action': 'corrected',
            'stored_inchikey': stored_ik,
            'pka_inchikey': result['inchikey'],
            'stored_smiles': stored_smiles,
            'pka_smiles': result['smiles'],
            'stored_formula': stored.get('formula', ''),
            'pka_formula': formula or '',
            'stored_charge': str(stored_charge),
            'pka_charge': str(target_charge), 'borderline_pka': ''})
        stats['pka_corrected'] += 1

    if correction_rows:
        fieldnames = ['cpd_id', 'compound_name', 'action',
                      'stored_inchikey', 'pka_inchikey',
                      'stored_smiles', 'pka_smiles',
                      'stored_formula', 'pka_formula',
                      'stored_charge', 'pka_charge', 'borderline_pka']
        with open(PROTONATION_CORRECTIONS_FILE, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in sorted(correction_rows, key=lambda r: r['cpd_id']):
                writer.writerow(row)

    logger.info("  pKa already correct: %d", stats['already_correct'])
    for key in ('pka_corrected', 'borderline_skipped', 'no_pka_data',
                'skipped_crossval'):
        if stats[key]:
            logger.info("  pKa %s: %d", key.replace('_', ' '), stats[key])
    return corrections


# ═══════════════════════════════════════════════════════════════════════════════
# Section 8: Apply Corrections & Write Dual-Format Output
# ═══════════════════════════════════════════════════════════════════════════════

def apply_corrections_dual_format(corrections, consistency_fixes, structures,
                                  all_file_path=None, unique_file_path=None,
                                  corrections_log_path=None):
    if all_file_path is None:
        all_file_path = ALL_STRUCTURES_FILE
    if unique_file_path is None:
        unique_file_path = UNIQUE_STRUCTURES_FILE
    if corrections_log_path is None:
        corrections_log_path = CORRECTIONS_LOG

    logger.info("Phase 6: Apply corrections & write output")
    logger.info("  All file: %s", all_file_path)
    logger.info("  Unique file: %s", unique_file_path)

    # Merge consistency fixes into corrections
    all_corrections = dict(corrections)
    for cpd_id, fixes in consistency_fixes.items():
        if cpd_id not in all_corrections:
            corr = {}
            stored = structures.get(cpd_id, {})
            for field, (_, new_val) in fixes.items():
                corr[field] = new_val
            if corr:
                corr.setdefault('smiles', stored.get('smiles', ''))
                corr.setdefault('inchi', stored.get('inchi', ''))
                corr.setdefault('inchikey', stored.get('inchikey', ''))
                corr['pubchem_cid'] = 'consistency'
                corr['result_type'] = 'CONSISTENCY_FIX'
                corr['strategy'] = 'phase0_consistency'
                corr['query'] = 'internal_validation'
                all_corrections[cpd_id] = corr

    # Read ALL rows into memory
    all_rows = []
    charged_indices = defaultdict(list)
    with open(all_file_path, "r") as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip('\n').split('\t')
            all_rows.append(cols)
            if len(cols) >= 8 and cols[2] == 'Charged':
                charged_indices[(cols[0], cols[1])].append(i)
    logger.info("  All file rows: %d", len(all_rows))

    # Compute formula/charge updates
    formula_updates = {}
    for cpd_id, corr in all_corrections.items():
        new_inchi = corr.get("inchi", "")
        if new_inchi:
            new_formula, new_charge = compute_formula_charge_from_inchi(
                new_inchi)
            if new_formula is not None:
                formula_updates[cpd_id] = (
                    new_formula,
                    str(new_charge) if new_charge is not None else "0")

    # Apply corrections and log changes
    log_header = ["timestamp", "cpd_id", "result_type", "field",
                  "old_value", "new_value", "pubchem_cid", "strategy",
                  "query"]
    ts = datetime.now().isoformat()
    total_changes = 0
    corrected_cpds = set()
    type_to_field = {"SMILE": "smiles", "InChIKey": "inchikey",
                     "InChI": "inchi"}

    with open(corrections_log_path, "w", newline="") as log_fh:
        log_writer = csv.writer(log_fh, delimiter="\t")
        log_writer.writerow(log_header)

        for (cpd_id, typ), indices in charged_indices.items():
            if cpd_id not in all_corrections:
                continue
            corr = all_corrections[cpd_id]
            field_key = type_to_field.get(typ)
            if not field_key:
                continue
            new_val = corr.get(field_key, "")
            if not new_val:
                continue
            for idx in indices:
                row = all_rows[idx]
                old_val = row[7]
                if old_val != new_val:
                    row[7] = new_val
                    if cpd_id in formula_updates:
                        row[5], row[6] = formula_updates[cpd_id]
                    log_writer.writerow([
                        ts, cpd_id, corr.get("result_type", ""), typ,
                        old_val, new_val, corr.get("pubchem_cid", ""),
                        corr.get("strategy", ""), corr.get("query", "")])
                    total_changes += 1
                    corrected_cpds.add(cpd_id)

        for cpd_id, (new_formula, new_charge) in formula_updates.items():
            old = structures.get(cpd_id, {})
            if old.get("formula", "") != new_formula or \
                    old.get("charge", "") != new_charge:
                # Apply formula/charge to all_rows
                for typ in ("SMILE", "InChI", "InChIKey"):
                    key = (cpd_id, typ)
                    if key in charged_indices:
                        for idx in charged_indices[key]:
                            all_rows[idx][5] = new_formula
                            all_rows[idx][6] = new_charge
                corr = all_corrections.get(cpd_id, {})
                log_writer.writerow([
                    ts, cpd_id, corr.get("result_type", ""),
                    "Formula/Charge",
                    f"{old.get('formula', '')}/{old.get('charge', '')}",
                    f"{new_formula}/{new_charge}",
                    corr.get("pubchem_cid", ""),
                    corr.get("strategy", ""), corr.get("query", "")])
                total_changes += 1
                corrected_cpds.add(cpd_id)

    # Canonicalize SMILES
    smiles_canonicalized = 0
    for (cpd_id, typ), indices in charged_indices.items():
        if typ != "SMILE":
            continue
        for idx in indices:
            row = all_rows[idx]
            old_smi = row[7]
            if not old_smi or old_smi == "null" or '*' in old_smi:
                continue
            try:
                mol = Chem.MolFromSmiles(old_smi)
                if mol:
                    canon = Chem.MolToSmiles(mol)
                    if canon and canon != old_smi:
                        row[7] = canon
                        smiles_canonicalized += 1
            except Exception:
                pass
    logger.info("  SMILES canonicalized: %d", smiles_canonicalized)

    # Write All file (atomic)
    tmp_all = all_file_path + ".tmp"
    with open(tmp_all, 'w') as fh:
        for row in all_rows:
            fh.write('\t'.join(row) + '\n')
    os.replace(tmp_all, all_file_path)
    logger.info("  All file written: %s", all_file_path)

    # Build and write Unique file
    unique = {}
    cpd_order = []
    seen_cpds = set()
    for row in all_rows:
        if len(row) < 8 or row[2] != "Charged":
            continue
        cpd_id, typ = row[0], row[1]
        alias_id, formula, charge, structure = row[3], row[5], row[6], row[7]
        key = (cpd_id, typ)
        if key not in unique:
            unique[key] = {
                "aliases": set(),
                "formula": formula if formula != "null" else "",
                "charge": charge if charge != "null" else "",
                "structure": structure if structure != "null" else ""}
        unique[key]["aliases"].add(alias_id)
        if not unique[key]["formula"] and formula != "null":
            unique[key]["formula"] = formula
            unique[key]["charge"] = charge if charge != "null" else ""
        if cpd_id not in seen_cpds:
            cpd_order.append(cpd_id)
            seen_cpds.add(cpd_id)

    # Propagate formula/charge across types
    for cpd_id in cpd_order:
        best_formula = best_charge = ""
        for typ in ("SMILE", "InChI", "InChIKey"):
            key = (cpd_id, typ)
            if key in unique and unique[key]["formula"]:
                best_formula = unique[key]["formula"]
                best_charge = unique[key]["charge"]
                break
        for typ in ("SMILE", "InChI", "InChIKey"):
            key = (cpd_id, typ)
            if key in unique and not unique[key]["formula"]:
                unique[key]["formula"] = best_formula
                unique[key]["charge"] = best_charge

    tmp_unique = unique_file_path + ".tmp"
    with open(tmp_unique, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["ID", "Type", "Aliases", "Formula", "Charge",
                         "Structure"])
        for cpd_id in cpd_order:
            for typ in ("SMILE", "InChIKey", "InChI"):
                key = (cpd_id, typ)
                if key not in unique:
                    continue
                entry = unique[key]
                writer.writerow([
                    cpd_id, typ, ";".join(sorted(entry["aliases"])),
                    entry["formula"], entry["charge"], entry["structure"]])
    os.replace(tmp_unique, unique_file_path)
    logger.info("  Unique file written: %s", unique_file_path)
    logger.info("  Compounds corrected: %d, fields changed: %d",
                len(corrected_cpds), total_changes)
    logger.info("  Corrections log: %s", corrections_log_path)

    if corrected_cpds:
        _run_post_correction_checks(corrected_cpds, unique)


def _run_post_correction_checks(corrected_cpds, unique):
    logger.info("  Running post-correction consistency checks...")
    inconsistencies = 0
    for cpd_id in sorted(corrected_cpds):
        cpd_data = {typ: unique[(cpd_id, typ)]["structure"]
                    for typ in ("SMILE", "InChI", "InChIKey")
                    if (cpd_id, typ) in unique}
        smiles = cpd_data.get("SMILE", "")
        inchi = cpd_data.get("InChI", "")
        inchikey = cpd_data.get("InChIKey", "")
        issues = []
        smi_mol = Chem.MolFromSmiles(smiles) if smiles else None
        if smiles and smi_mol is None:
            issues.append("SMILES fails RDKit parsing")
        if inchi and MolFromInchi(inchi) is None:
            issues.append("InChI fails RDKit parsing")
        if inchi and inchikey:
            try:
                computed = InchiToInchiKey(inchi)
                if computed and computed != inchikey:
                    issues.append(f"InChIKey mismatch: stored="
                                  f"{inchikey[:14]} computed={computed[:14]}")
            except Exception:
                issues.append("InChIKey computation failed")
        if smi_mol and inchi:
            try:
                key_smi = Chem.inchi.MolToInchiKey(smi_mol)
                key_inchi = InchiToInchiKey(inchi)
                if (key_smi and key_inchi
                        and key_smi.split("-")[0] != key_inchi.split("-")[0]):
                    issues.append("SMILES/InChI connectivity mismatch")
            except Exception:
                pass
        if issues:
            inconsistencies += 1
            logger.warning("  Consistency issue %s: %s",
                           cpd_id, "; ".join(issues))
    if inconsistencies:
        logger.warning("  %d corrected compounds have consistency issues",
                       inconsistencies)
    else:
        logger.info("  All %d corrected compounds pass consistency checks",
                    len(corrected_cpds))


# ═══════════════════════════════════════════════════════════════════════════════
# Section 9: Reporting
# ═══════════════════════════════════════════════════════════════════════════════

def generate_report(conn, db_lock, candidates, structures, names,
                    report_file, mismatch_file, protonation_file, stereo_file,
                    tautomer_file=None):
    logger.info("Phase 4: Classification & reporting")
    cur = conn.execute("SELECT * FROM cache WHERE cpd_id IN ({})".format(
        ",".join("?" * len(candidates))), candidates)

    report_header = [
        "cpd_id", "compound_name", "name_queried", "strategy", "result_type",
        "stored_inchikey", "pubchem_inchikey", "pubchem_cid",
        "stored_smiles", "pubchem_smiles"]
    mismatch_header = report_header + ["mismatch_detail"]
    stereo_header = report_header + [
        "stored_stereo_defined", "stored_stereo_potential",
        "pubchem_stereo_defined", "pubchem_stereo_potential"]

    stats = defaultdict(int)
    mismatch_substats = defaultdict(int)
    report_rows, mismatch_work = [], []
    protonation_rows, stereo_rows = [], []
    reclassifications = {}

    status_to_rt = {"not_found": "NOT_FOUND", "ambiguous": "AMBIGUOUS",
                    "xref_conflict": "XREF_CONFLICT"}

    for row in cur.fetchall():
        cpd_id = row[0]
        strategy, query, status = row[1], row[2], row[3]
        pub_cid, pub_ik, pub_smiles = row[4] or "", row[5] or "", row[6] or ""
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        stored_smiles = structures.get(cpd_id, {}).get("smiles", "")
        cpd_names = names.get(cpd_id, []) if names else []
        compound_name = cpd_names[0] if cpd_names else ""

        rt = status_to_rt.get(status) or compare_inchikeys(stored_ik, pub_ik)
        stats[rt] += 1
        out_row = [cpd_id, compound_name, query, strategy, rt,
                   stored_ik, pub_ik, pub_cid, stored_smiles, pub_smiles]
        report_rows.append(out_row)

        if rt == "MISMATCH":
            stored_struct = structures.get(cpd_id, {})
            pub_formula = ""
            if pub_smiles:
                try:
                    p_mol = Chem.MolFromSmiles(pub_smiles)
                    if p_mol:
                        pub_formula = rdMolDescriptors.CalcMolFormula(p_mol)
                except Exception:
                    pass
            mismatch_work.append((cpd_id, stored_struct,
                                  {"smiles": pub_smiles, "inchi": "",
                                   "inchikey": pub_ik,
                                   "formula": pub_formula}, out_row))
        elif rt == "PROTONATION_DIFF":
            protonation_rows.append(out_row)
        elif rt == "STEREO_DIFF":
            s_spec, s_pot = count_defined_stereo(stored_smiles)
            p_spec, p_pot = count_defined_stereo(pub_smiles)
            stereo_rows.append(out_row + [
                s_spec if s_spec is not None else "",
                s_pot if s_pot is not None else "",
                p_spec if p_spec is not None else "",
                p_pot if p_pot is not None else ""])

    # Classify mismatches in parallel
    mismatch_rows = []
    if mismatch_work:
        num_workers = min(os.cpu_count() or 4, 32)
        logger.info("  Classifying %d mismatches...", len(mismatch_work))
        with Pool(num_workers) as pool:
            for cpd_id, out_row, detail, reclassify in pool.imap_unordered(
                    classify_mismatch_worker, mismatch_work, chunksize=32):
                subcat = detail.split("(")[0].strip() if "(" in detail \
                    else detail
                mismatch_substats[subcat] += 1
                mismatch_rows.append(out_row + [detail])
                if reclassify:
                    reclassifications[cpd_id] = reclassify

    # Reclassify eligible MISMATCH sub-types in cache
    if reclassifications:
        cpd_ids_list = list(reclassifications.keys())
        with db_lock:
            cur_r = conn.execute(
                "SELECT cpd_id, strategy, query, pubchem_cid, pubchem_inchikey, "
                "pubchem_smiles FROM cache WHERE cpd_id IN ({})".format(
                    ",".join("?" * len(cpd_ids_list))), cpd_ids_list)
            cache_rows = {r[0]: r for r in cur_r.fetchall()}
        reclassify_rows = []
        for cpd_id, new_rt in reclassifications.items():
            r = cache_rows.get(cpd_id)
            if r:
                reclassify_rows.append((cpd_id, r[1], r[2], "found",
                                        r[3], r[4], r[5], time.time()))
            stats["MISMATCH"] -= 1
            stats[new_rt] += 1
            for out_row in report_rows:
                if out_row[0] == cpd_id:
                    out_row[4] = new_rt
                    pub_smiles = out_row[9]
                    stored_smiles = out_row[8]
                    if new_rt == "PROTONATION_DIFF":
                        protonation_rows.append(out_row)
                    elif new_rt == "STEREO_DIFF":
                        s_spec, s_pot = count_defined_stereo(stored_smiles)
                        p_spec, p_pot = count_defined_stereo(pub_smiles)
                        stereo_rows.append(out_row + [
                            s_spec if s_spec is not None else "",
                            s_pot if s_pot is not None else "",
                            p_spec if p_spec is not None else "",
                            p_pot if p_pot is not None else ""])
                    break
        if reclassify_rows:
            save_batch_to_cache(conn, db_lock, reclassify_rows)
        stereo_n = sum(1 for v in reclassifications.values()
                       if v == "STEREO_DIFF")
        prot_n = sum(1 for v in reclassifications.values()
                     if v == "PROTONATION_DIFF")
        logger.info("  Reclassified: %d (STEREO: %d, PROTONATION: %d)",
                    len(reclassifications), stereo_n, prot_n)

    # Write reports
    _write_tsv(report_file, report_header, report_rows)
    _write_tsv(mismatch_file, mismatch_header, mismatch_rows)
    if protonation_rows:
        _write_tsv(protonation_file, report_header, protonation_rows)
    if stereo_rows:
        _write_tsv(stereo_file, stereo_header, stereo_rows)
    if tautomer_file and mismatch_rows:
        taut_rows = [r for r in mismatch_rows if 'tautomer' in r[-1].lower()]
        if taut_rows:
            _write_tsv(tautomer_file, mismatch_header, taut_rows)

    logger.info("=" * 50)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 50)
    total = sum(stats.values())
    logger.info("Total compounds checked: %d", total)
    for key in ("MATCH", "STEREO_DIFF", "PROTONATION_DIFF", "MISMATCH",
                "NOT_FOUND", "NO_KEY", "AMBIGUOUS", "XREF_CONFLICT"):
        if stats[key]:
            pct = 100.0 * stats[key] / total if total else 0
            logger.info("  %-20s: %6d (%5.1f%%)", key, stats[key], pct)
    if mismatch_substats:
        logger.info("  MISMATCH sub-categories:")
        for subcat, cnt in sorted(mismatch_substats.items(),
                                  key=lambda x: -x[1]):
            logger.info("    %-40s: %4d", subcat, cnt)
    return reclassifications


def generate_comparison_images(conn, candidates, structures, images_dir,
                               max_per_type=20, accepted_cpds=None):
    from PIL import Image, ImageDraw, ImageFont
    from rdkit.Chem import Draw

    os.makedirs(images_dir, exist_ok=True)
    if accepted_cpds is None:
        accepted_cpds = set()
    for f in os.listdir(images_dir):
        if f.endswith(".png"):
            os.remove(os.path.join(images_dir, f))

    cur = conn.execute(
        "SELECT cpd_id, strategy, query, status, pubchem_cid, "
        "pubchem_inchikey FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))), candidates)
    by_type = {"MISMATCH": [], "STEREO_DIFF": [], "PROTONATION_DIFF": []}
    cid_map = {}
    for row in cur.fetchall():
        cpd_id, _, _, status, pub_cid, pub_ik = row
        if status != "found" or not pub_cid:
            continue
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        rt = compare_inchikeys(stored_ik, pub_ik)
        if rt in by_type:
            by_type[rt].append(cpd_id)
            cid_map[cpd_id] = pub_cid

    # Batch-fetch all needed SMILES
    pubchem_smiles = {}
    all_cpds = [c for cpds in by_type.values() for c in cpds]
    if all_cpds:
        cur2 = conn.execute(
            "SELECT cpd_id, pubchem_smiles FROM corrections "
            "WHERE cpd_id IN ({})".format(",".join("?" * len(all_cpds))),
            all_cpds)
        pubchem_smiles = {r[0]: r[1] for r in cur2.fetchall()}
    missing_cids = {}
    for c in all_cpds:
        if c not in pubchem_smiles and c in cid_map:
            missing_cids.setdefault(cid_map[c], []).append(c)
    if missing_cids:
        cid_list = sorted(missing_cids.keys())
        for i in range(0, len(cid_list), CID_BATCH_SIZE):
            batch = cid_list[i:i + CID_BATCH_SIZE]
            for cid, props in query_cid_properties_batch(batch).items():
                if props.get("smiles"):
                    for cpd_id in missing_cids.get(cid, []):
                        pubchem_smiles[cpd_id] = props["smiles"]

    try:
        font_large = ImageFont.truetype("DejaVuSansMono.ttf", 16)
        font_small = ImageFont.truetype("DejaVuSansMono.ttf", 12)
    except (OSError, IOError):
        try:
            font_large = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 16)
            font_small = ImageFont.truetype(
                "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", 12)
        except (OSError, IOError):
            font_large = font_small = ImageFont.load_default()

    def _render_image(cpd_id, result_type, tag, tag_color, suffix):
        stored_smi = structures.get(cpd_id, {}).get("smiles", "")
        pub_smi = pubchem_smiles.get(cpd_id, "")
        mol_ms = Chem.MolFromSmiles(stored_smi) if stored_smi else None
        mol_pc = Chem.MolFromSmiles(pub_smi) if pub_smi else None
        if mol_ms is None and mol_pc is None:
            return False
        mol_size = (350, 300)
        img_ms = (Draw.MolToImage(mol_ms, size=mol_size) if mol_ms
                  else Image.new("RGB", mol_size, "white"))
        img_pc = (Draw.MolToImage(mol_pc, size=mol_size) if mol_pc
                  else Image.new("RGB", mol_size, "white"))
        margin, header_h, gap = 10, 50, 20
        w = mol_size[0] * 2 + gap + margin * 2
        h = mol_size[1] + header_h + 45 + margin * 2
        canvas = Image.new("RGB", (w, h), "white")
        draw = ImageDraw.Draw(canvas)
        draw.text((margin, 5), f"{cpd_id}  [{result_type}]",
                  fill="black", font=font_large)
        draw.text((w - margin - len(tag) * 10, 5), tag,
                  fill=tag_color, font=font_large)
        left_x, right_x = margin, margin + mol_size[0] + gap
        draw.text((left_x + mol_size[0] // 2 - 40, header_h - 18),
                  "ModelSEED", fill="blue", font=font_large)
        draw.text((right_x + mol_size[0] // 2 - 35, header_h - 18),
                  "PubChem", fill="red", font=font_large)
        canvas.paste(img_ms, (left_x, header_h))
        canvas.paste(img_pc, (right_x, header_h))
        prefix = {"MISMATCH": "mismatch", "STEREO_DIFF": "stereo",
                  "PROTONATION_DIFF": "prot"}[result_type]
        canvas.save(os.path.join(images_dir,
                                 f"{prefix}_{suffix}_{cpd_id}.png"))
        return True

    total_images = 0
    for result_type, cpd_list in by_type.items():
        accepted = sorted(c for c in cpd_list if c in accepted_cpds)
        not_accepted = sorted(c for c in cpd_list if c not in accepted_cpds)
        rej_tag = ("NOT_AUTO_CORRECTED" if result_type == "PROTONATION_DIFF"
                   else "REJECTED")
        count = 0
        for cpd_id in accepted:
            if count >= max_per_type:
                break
            if _render_image(cpd_id, result_type, "CORRECTED", "green",
                             "accepted"):
                count += 1
                total_images += 1
        for cpd_id in not_accepted:
            if count >= max_per_type:
                break
            if _render_image(cpd_id, result_type, rej_tag, "red",
                             "rejected"):
                count += 1
                total_images += 1
    logger.info("  Structure images: %d (in %s/)", total_images, images_dir)


def generate_different_compound_review(conn, candidates, structures, names,
                                       mismatch_file, review_file, workers=5):
    review_cpds = {}
    if not os.path.exists(mismatch_file):
        return
    with open(mismatch_file, 'r') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if 'REVIEW:different_compound' in row.get('mismatch_detail', ''):
                review_cpds[row['cpd_id']] = row
    if not review_cpds:
        return
    logger.info("Generating review report for %d different-compound "
                "mismatches", len(review_cpds))

    to_lookup = {structures.get(c, {}).get('inchikey', ''): c
                 for c in review_cpds
                 if structures.get(c, {}).get('inchikey', '')}
    correct_cids = {}
    if to_lookup:
        ik_to_cpds = defaultdict(list)
        for ik, cpd_id in to_lookup.items():
            ik_to_cpds[ik].append(cpd_id)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(query_inchikey, ik): ik
                       for ik in ik_to_cpds}
            for future in as_completed(futures):
                ik = futures[future]
                result = future.result()
                if result is not None:
                    cid, pub_ik, smiles = result
                    for cpd_id in ik_to_cpds[ik]:
                        correct_cids[cpd_id] = (cid, pub_ik, smiles)

    review_rows = []
    for cpd_id, mr in review_cpds.items():
        stored = structures.get(cpd_id, {})
        strategy = mr.get('strategy', '')
        bad_cid = mr.get('pubchem_cid', '')
        if 'kegg_xref' in strategy:
            cause = 'wrong_kegg_mapping'
        elif 'chebi_xref' in strategy:
            cause = 'wrong_chebi_mapping'
        elif 'name_lookup' in strategy:
            cause = 'ambiguous_name'
        else:
            cause = 'unknown'

        correct = correct_cids.get(cpd_id)
        if correct:
            c_cid, c_ik = correct[0], correct[1]
            action = (f"Same CID {bad_cid}" if c_cid == bad_cid
                      else f"Correct CID is {c_cid}; xref points to "
                           f"wrong CID {bad_cid}")
        else:
            c_cid, c_ik = '', ''
            action = (f"Not found in PubChem; {strategy} xref maps to "
                      f"wrong compound")

        detail = mr.get('mismatch_detail', '')
        tanimoto = ''
        if 'tanimoto=' in detail:
            try:
                tanimoto = detail.split('tanimoto=')[1].split(',')[0] \
                    .split(')')[0]
            except (IndexError, ValueError):
                pass

        cpd_names = names.get(cpd_id, []) if names else []
        review_rows.append({
            'cpd_id': cpd_id,
            'compound_name': cpd_names[0] if cpd_names else "",
            'strategy': strategy, 'query': mr.get('name_queried', ''),
            'likely_cause': cause, 'suggested_action': action,
            'stored_inchikey': stored.get('inchikey', ''),
            'stored_smiles': stored.get('smiles', ''),
            'stored_formula': stored.get('formula', ''),
            'pubchem_cid_bad': bad_cid,
            'pubchem_inchikey_bad': mr.get('pubchem_inchikey', ''),
            'pubchem_smiles_bad': mr.get('pubchem_smiles', ''),
            'correct_pubchem_cid': c_cid,
            'correct_pubchem_inchikey': c_ik,
            'tanimoto_similarity': tanimoto})

    review_rows.sort(
        key=lambda r: (0 if r['correct_pubchem_cid'] else 1, r['cpd_id']))
    fieldnames = ['cpd_id', 'compound_name', 'strategy', 'query',
                  'likely_cause', 'suggested_action',
                  'stored_inchikey', 'stored_smiles', 'stored_formula',
                  'pubchem_cid_bad', 'pubchem_inchikey_bad',
                  'pubchem_smiles_bad', 'correct_pubchem_cid',
                  'correct_pubchem_inchikey', 'tanimoto_similarity']
    with open(review_file, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(review_rows)
    logger.info("  Review report: %s (%d entries)",
                review_file, len(review_rows))


# ═══════════════════════════════════════════════════════════════════════════════
# Section 10: Tests & Main
# ═══════════════════════════════════════════════════════════════════════════════

def run_tests():
    print("Running unit tests...")
    passed = failed = 0
    tests = [
        ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "XLYOFNOQVPJJNP-UHFFFAOYSA-N",
         "MATCH"),
        ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "XLYOFNOQVPJJNP-UHFFFAOYSA-M",
         "PROTONATION_DIFF"),
        ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "XLYOFNOQVPJJNP-ABCDEFGHIJ-N",
         "STEREO_DIFF"),
        ("XLYOFNOQVPJJNP-UHFFFAOYSA-N", "ZZZZZZZZZZZZZ-UHFFFAOYSA-N",
         "MISMATCH"),
        ("", "XLYOFNOQVPJJNP-UHFFFAOYSA-N", "NO_KEY"),
    ]
    for stored, pubchem, expected in tests:
        result = compare_inchikeys(stored, pubchem)
        if result == expected:
            passed += 1
        else:
            print(f"  FAIL: compare_inchikeys({stored[:14]}, "
                  f"{pubchem[:14]}) = {result}, expected {expected}")
            failed += 1

    f, c = compute_formula_charge_from_inchi("InChI=1S/H2O/h1H2")
    if f == "H2O" and c == 0:
        passed += 1
    else:
        print(f"  FAIL: compute_formula_charge_from_inchi(water) = ({f}, {c})")
        failed += 1

    mol = Chem.MolFromSmiles("C1=CC=CC=C1")
    if mol and Chem.MolToSmiles(mol) == "c1ccccc1":
        passed += 1
    else:
        print("  FAIL: canonical benzene")
        failed += 1

    d, p, net = predict_charge_from_pka(
        [(1, 1, 3.0), (1, 2, 10.0)], [(1, 3, 9.0)])
    if d == 1 and p == 1 and net == 0:
        passed += 1
    else:
        print(f"  FAIL: predict_charge_from_pka = ({d}, {p}, {net})")
        failed += 1

    result = _classify_mismatch(
        {"smiles": "O", "inchi": "InChI=1S/H2O/h1H2",
         "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N", "formula": "H2O"},
        {"smiles": "O", "inchi": "InChI=1S/H2O/h1H2",
         "inchikey": "XLYOFNOQVPJJNP-UHFFFAOYSA-N", "formula": "H2O"})
    if isinstance(result, tuple) and len(result) == 2:
        passed += 1
    else:
        print(f"  FAIL: _classify_mismatch return type = {type(result)}")
        failed += 1

    print(f"\nTests: {passed} passed, {failed} failed")
    return failed == 0


def main():
    args = parse_args()

    if args.test:
        sys.exit(0 if run_tests() else 1)

    _file_handler = logging.FileHandler(LOG_FILE, mode="w")
    _file_handler.setLevel(logging.DEBUG)
    _file_handler.setFormatter(_log_formatter)
    logger.addHandler(_file_handler)
    _console_handler = logging.StreamHandler()
    _console_handler.setLevel(logging.INFO)
    _console_handler.setFormatter(_log_formatter)
    logger.addHandler(_console_handler)

    if args.clear_cache and os.path.exists(CACHE_DB):
        os.remove(CACHE_DB)
        logger.info("Cache cleared.")

    logger.info("Loading local data...")
    structures = load_structures(ALL_STRUCTURES_FILE)
    names = load_names(NAMES_FILE)
    external_ids = load_external_ids(ALIASES_FILE)

    total_cpds = len(structures)
    cpds_with_ik = sum(1 for s in structures.values() if 'inchikey' in s)
    cpds_with_names = sum(1 for c in structures if c in names)
    logger.info("Coverage statistics:")
    logger.info("  Total compounds: %d", total_cpds)
    logger.info("  With InChIKey: %d (%d without)",
                cpds_with_ik, total_cpds - cpds_with_ik)
    logger.info("  With names: %d (%d without)",
                cpds_with_names, total_cpds - cpds_with_names)

    consistency_fixes = run_phase0_consistency(
        structures, report_file=CONSISTENCY_FILE)

    pka_data = load_pka_from_db()
    logger.info("ChemAxon pKa data loaded: %d compounds", len(pka_data))

    candidates = [cpd_id for cpd_id in sorted(structures.keys())
                  if cpd_id in names and "inchikey" in structures[cpd_id]]
    logger.info("Validation candidates (name + InChIKey): %d", len(candidates))

    if args.limit > 0:
        candidates = candidates[:args.limit]
        logger.info("Limited to first %d compounds", args.limit)

    conn = init_cache(CACHE_DB)
    try:
        db_lock = threading.Lock()
        cached = get_cached(conn)
        to_process = [c for c in candidates if c not in cached]
        logger.info("Already cached: %d", len(candidates) - len(to_process))
        logger.info("To process: %d", len(to_process))

        if to_process:
            needs_name_lookup = run_phase1_pubchem_lookup(
                to_process, external_ids, structures, conn, db_lock,
                args.workers, candidates)
            run_phase2_name_lookup(needs_name_lookup, names, structures, conn,
                                   db_lock, args.workers)
            run_phase3_recovery(candidates, structures, conn, db_lock,
                                args.workers)
        else:
            _resolve_cached_xref_conflicts(
                conn, db_lock, candidates, structures, args.workers)

        reclassifications = generate_report(
            conn, db_lock, candidates, structures, names=names,
            report_file=REPORT_FILE, mismatch_file=MISMATCH_FILE,
            protonation_file=PROTONATION_FILE, stereo_file=STEREO_FILE,
            tautomer_file=TAUTOMER_FILE)

        corrections = {}
        if args.apply or args.images:
            corrections = run_phase5_corrections(
                conn, db_lock, candidates, structures, pka_data, names,
                args.workers, skip_ph7=args.skip_ph7,
                reclassifications=reclassifications)

        if args.apply:
            apply_corrections_dual_format(
                corrections, consistency_fixes, structures,
                all_file_path=ALL_STRUCTURES_FILE,
                unique_file_path=UNIQUE_STRUCTURES_FILE,
                corrections_log_path=CORRECTIONS_LOG)

        if args.images:
            generate_comparison_images(conn, candidates, structures,
                                       images_dir=IMAGES_DIR,
                                       max_per_type=args.max_images,
                                       accepted_cpds=set(corrections.keys()))

        if args.review_report:
            generate_different_compound_review(
                conn, candidates, structures, names,
                mismatch_file=MISMATCH_FILE, review_file=REVIEW_FILE,
                workers=args.workers)
    finally:
        conn.close()
    logger.info("Done.")


if __name__ == "__main__":
    main()
