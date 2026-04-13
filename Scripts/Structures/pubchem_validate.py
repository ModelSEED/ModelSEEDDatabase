#!/usr/bin/env python3
# usage: python pubchem_validate.py --apply --output Unique_ModelSEED_Structures_modified.txt
"""PubChem name-structure cross-validation for ModelSEED compounds.

Validates whether compound names match their stored structures by querying
PubChem for InChIKeys via ChEBI xref, KEGG xref, name lookup, and direct
InChIKey lookup, then comparing against locally stored InChIKeys.

With --apply, fetches full structures from PubChem for STEREO_DIFF compounds
and applies corrections only when they demonstrably improve the data:
  - MISMATCH: never auto-corrected (classified for manual review)
  - STEREO_DIFF: accepted only if PubChem has >= defined stereocenters
  - PROTONATION_DIFF: never auto-corrected (pH-dependent)

Safety features:
  - Timestamped backup of structures file before any modification
  - Corrections log is appended (not overwritten) with timestamps
  - Formula and Charge columns updated to match new structures
  - Cross-validation of xref sources (ChEBI vs KEGG)
  - Multi-name cross-validation (ambiguity detection)

Usage:
    python pubchem_validate.py              # full run (resumable)
    python pubchem_validate.py --limit 100  # test with first 100 compounds
    python pubchem_validate.py --apply      # full run + apply corrections
    python pubchem_validate.py --resume     # continue from cache (default)
    python pubchem_validate.py --workers 3  # adjust concurrency (default 10)
    python pubchem_validate.py --clear-cache  # clear cache for fresh run
"""

import argparse
import logging
import os
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from pubchem_api import (query_xref, query_inchikey, query_inchi,
                         query_cid_properties, query_cid_properties_batch,
                         CID_BATCH_SIZE,
                         _query_names_batch_recursive, NAME_BATCH_SIZE)
from data_io import (load_structures, load_names, load_external_ids,
                     init_cache, get_cached, save_batch_to_cache)
from structure_compare import compare_inchikeys, pick_best_names
from corrections import (fetch_corrections, apply_corrections,
                        normalize_corrections_ph7)
from reporting import (generate_report, generate_comparison_images,
                       generate_different_compound_review)
from protonation import run_phase4_protonation

# ---------------------------------------------------------------------------
# File path constants
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)

STRUCTURES_FILE = os.path.join(SCRIPT_DIR, "All_ModelSEED_Structures.txt")
NAMES_FILE = os.path.join(
    BASE_DIR, "ModelSEEDDatabase", "Biochemistry", "Aliases",
    "Unique_ModelSEED_Compound_Names.txt",
)
ALIASES_FILE = os.path.join(
    BASE_DIR, "ModelSEEDDatabase", "Biochemistry", "Aliases",
    "Unique_ModelSEED_Compound_Aliases.txt",
)
CACHE_DB = os.path.join(SCRIPT_DIR, "pubchem_cache_all.sqlite")
REPORT_FILE = os.path.join(SCRIPT_DIR, "pubchem_validation_report.tsv")
MISMATCH_FILE = os.path.join(SCRIPT_DIR, "pubchem_mismatches.tsv")
CORRECTIONS_LOG = os.path.join(SCRIPT_DIR, "pubchem_corrections_log.tsv")
PROTONATION_FILE = os.path.join(SCRIPT_DIR, "pubchem_protonation_diffs.tsv")
STEREO_FILE = os.path.join(SCRIPT_DIR, "pubchem_stereo_diffs.tsv")
REJECTED_FILE = os.path.join(SCRIPT_DIR, "pubchem_rejected_corrections.tsv")
PROTONATION_CORRECTIONS_FILE = os.path.join(
    SCRIPT_DIR, "pubchem_protonation_corrections.tsv")
REVIEW_FILE = os.path.join(
    SCRIPT_DIR, "pubchem_review_different_compounds.tsv")

LOG_FILE = os.path.join(SCRIPT_DIR, "pubchem_validate_run.log")
IMAGES_DIR = os.path.join(SCRIPT_DIR, "struct_imgs")

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
logger = logging.getLogger("pubchem_validate")
logger.setLevel(logging.DEBUG)
_log_formatter = logging.Formatter(
    "%(asctime)s  %(levelname)-8s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
_file_handler = logging.FileHandler(LOG_FILE, mode="w")
_file_handler.setLevel(logging.DEBUG)
_file_handler.setFormatter(_log_formatter)
_file_handler.stream = open(LOG_FILE, "w", buffering=1)  # line-buffered
logger.addHandler(_file_handler)


# ---------------------------------------------------------------------------
# Phase 1.7: Resolve cached xref conflicts by fetching CID properties
# ---------------------------------------------------------------------------

def run_phase17_resolve_cached_conflicts(conn, db_lock, candidates,
                                         structures, workers=5):
    """Phase 1.7: Resolve xref conflicts from cache by fetching CID properties.

    For each xref_conflict row:
    1. Parse CIDs from query string (format: chebi_xref->CID{X};kegg_xref={Y}->CID{Z})
    2. Fetch InChIKey/SMILES for each CID via query_cid_properties()
    3. Compare each to stored InChIKey using compare_inchikeys()
    4. Pick best match, update cache row

    Returns (resolved_count, still_conflicting_count).
    """
    import re

    cur = conn.execute(
        "SELECT cpd_id, query FROM cache WHERE status='xref_conflict' "
        "AND cpd_id IN ({})".format(",".join("?" * len(candidates))),
        candidates,
    )
    conflicts = cur.fetchall()
    if not conflicts:
        return 0, 0

    logger.info("Phase 1.7: Resolving %d cached xref conflicts", len(conflicts))

    # Parse CIDs from query strings and deduplicate
    cpd_cids = {}  # cpd_id -> [(cid, source_hint), ...]
    all_cids = set()
    for cpd_id, query_str in conflicts:
        cids_found = re.findall(r'CID(\d+)', query_str)
        pairs = []
        for cid_str in cids_found:
            all_cids.add(cid_str)
            # Determine source hint from position in query string
            idx = query_str.index(f'CID{cid_str}')
            prefix = query_str[:idx]
            source = 'chebi' if 'chebi' in prefix.split(';')[-1] else 'kegg'
            pairs.append((cid_str, source))
        cpd_cids[cpd_id] = pairs

    logger.info("  Unique CIDs to fetch: %d", len(all_cids))

    # Fetch properties for all unique CIDs in batches
    cid_props = {}  # cid -> {smiles, inchi, inchikey} or None
    cid_list = sorted(all_cids)
    batches = [cid_list[i:i + CID_BATCH_SIZE]
               for i in range(0, len(cid_list), CID_BATCH_SIZE)]

    logger.info("  Batches: %d (batch size %d)", len(batches), CID_BATCH_SIZE)

    if tqdm:
        pbar = tqdm(total=len(batches), desc="Phase 1.7: CID batches")

    batch_done = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(query_cid_properties_batch, batch): batch
                   for batch in batches}
        for future in as_completed(futures):
            batch_result = future.result()
            cid_props.update(batch_result)
            batch_done += 1
            if tqdm:
                pbar.update(1)
            elif batch_done % 10 == 0:
                logger.info("  Phase 1.7 progress: %d/%d batches done",
                            batch_done, len(batches))

    if tqdm:
        pbar.close()

    logger.info("  CIDs fetched: %d", len(cid_props))

    # Resolve each conflict
    _result_priority = {
        "MATCH": 0, "PROTONATION_DIFF": 1, "STEREO_DIFF": 2,
        "MISMATCH": 3, "NO_KEY": 4,
    }
    resolution_rows = []
    resolved_count = 0
    still_conflicting = 0

    for cpd_id, query_str in conflicts:
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        if not stored_ik:
            still_conflicting += 1
            continue

        scored = []
        for cid_str, source in cpd_cids[cpd_id]:
            props = cid_props.get(cid_str)
            if props is None or not props.get("inchikey"):
                continue
            result_type = compare_inchikeys(stored_ik, props["inchikey"])
            priority = _result_priority.get(result_type, 99)
            source_pref = 0 if source == 'chebi' else 1
            scored.append((priority, source_pref, cid_str, props, result_type))

        if not scored:
            still_conflicting += 1
            continue

        scored.sort()
        best_priority, _, best_cid, best_props, best_rt = scored[0]

        if best_priority <= 2:  # MATCH, PROTONATION_DIFF, or STEREO_DIFF
            query_resolved = (f"xref_resolved:CID{best_cid}({best_rt})")
            resolution_rows.append((
                cpd_id, "xref_resolved", query_resolved, "found",
                best_cid, best_props["inchikey"], best_props.get("smiles", ""),
                time.time()))
            resolved_count += 1
        else:
            still_conflicting += 1

    if resolution_rows:
        save_batch_to_cache(conn, db_lock, resolution_rows)

    logger.info("  Resolved: %d", resolved_count)
    logger.info("  Still conflicting: %d", still_conflicting)
    return resolved_count, still_conflicting


# ---------------------------------------------------------------------------
# Phase 1: xref lookups (ChEBI + KEGG) with cross-validation
# ---------------------------------------------------------------------------

def run_phase1_xref(to_process, external_ids, structures, conn, db_lock,
                    workers):
    """Query PubChem by ChEBI and KEGG xrefs with cross-validation.

    Returns (xref_resolved_cpds, xref_conflicts, xref_conflict_pairs,
             needs_name_lookup).
    """
    chebi_to_cpds = defaultdict(list)
    kegg_to_cpds = defaultdict(list)
    cpd_has_chebi = set()
    cpd_has_kegg = set()
    needs_name_lookup = []

    for cpd_id in to_process:
        ext = external_ids.get(cpd_id, {"chebi": [], "kegg": []})
        has_chebi = False
        has_kegg = False
        for cid_val in ext["chebi"]:
            chebi_key = f"CHEBI:{cid_val}"
            chebi_to_cpds[chebi_key].append(cpd_id)
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
    both_count = len(cpd_has_chebi & cpd_has_kegg)

    logger.info("Phase 1: xref lookups")
    logger.info("  Unique ChEBI IDs: %d", len(unique_chebi))
    logger.info("  Unique KEGG IDs: %d", len(unique_kegg))
    logger.info("  Compounds with both: %d", both_count)
    logger.info("  Total xref queries (deduplicated): %d", total_xref_queries)
    logger.info("  Compounds needing name lookup: %d", len(needs_name_lookup))

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
        result = query_xref(xref_id)
        return (xref_id, strategy, result)

    if tqdm:
        pbar = tqdm(total=total_xref_queries, desc="Phase 1: xref queries")

    xref_done = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = []
        for cid in unique_chebi:
            futures.append(executor.submit(do_xref_query, cid, "chebi_xref"))
        for kid in unique_kegg:
            futures.append(executor.submit(do_xref_query, kid, "kegg_xref"))

        for future in as_completed(futures):
            xref_id, strategy, result = future.result()
            xref_done += 1
            if tqdm:
                pbar.update(1)
            if xref_done % 500 == 0:
                logger.info("  Phase 1 progress: %d/%d xref queries done",
                            xref_done, total_xref_queries)

            if result is None:
                continue

            cid, inchikey, smiles = result
            if strategy == "chebi_xref":
                cpd_ids = chebi_to_cpds[xref_id]
            else:
                cpd_ids = kegg_to_cpds[xref_id]

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
                                    row = conn.execute(
                                        "SELECT pubchem_cid, "
                                        "pubchem_inchikey, pubchem_smiles, "
                                        "strategy, query FROM cache "
                                        "WHERE cpd_id=?",
                                        (cpd_id,)).fetchone()
                                if row:
                                    xref_conflict_pairs[cpd_id].append(
                                        (row[0], row[1], row[2],
                                         row[3], row[4]))
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
    if tqdm:
        pbar.close()

    logger.info("  Resolved via xref: %d", len(xref_resolved_cpds))
    if xref_conflicts:
        logger.info("  xref conflicts (ChEBI/KEGG disagree): %d",
                    len(xref_conflicts))

    return (xref_resolved_cpds, xref_conflicts, xref_conflict_pairs,
            needs_name_lookup)


# ---------------------------------------------------------------------------
# Phase 1.5: Resolve xref conflicts using InChIKey comparison
# ---------------------------------------------------------------------------

def run_phase15_resolve_conflicts(xref_conflict_pairs, structures, conn,
                                  db_lock):
    """Resolve xref conflicts by comparing candidates to stored InChIKeys.

    Returns (resolved_cpds, still_conflicting_cpds) as sets.
    """
    if not xref_conflict_pairs:
        return set(), set()

    logger.info("Phase 1.5: Resolving xref conflicts")
    logger.info("  Conflicts with candidate data: %d",
                len(xref_conflict_pairs))
    resolved_count = 0
    still_conflicting = 0
    _result_priority = {
        "MATCH": 0, "PROTONATION_DIFF": 1, "STEREO_DIFF": 2,
        "MISMATCH": 3, "NO_KEY": 4,
    }
    resolution_rows = []
    resolved_cpds = set()
    still_conflicting_cpds = set()

    for cpd_id, candidates in xref_conflict_pairs.items():
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        if not stored_ik:
            still_conflicting += 1
            still_conflicting_cpds.add(cpd_id)
            continue

        scored = []
        for c_cid, c_ik, c_smi, c_strategy, c_xref_id in candidates:
            if not c_ik:
                continue
            result_type = compare_inchikeys(stored_ik, c_ik)
            priority = _result_priority.get(result_type, 99)
            source_pref = 0 if "chebi" in c_strategy else 1
            scored.append((priority, source_pref, c_cid, c_ik, c_smi,
                           c_strategy, c_xref_id, result_type))

        if not scored:
            still_conflicting += 1
            still_conflicting_cpds.add(cpd_id)
            continue

        scored.sort()
        best = scored[0]
        best_priority, _, b_cid, b_ik, b_smi, b_strat, b_xref, b_rt = best

        if best_priority <= 2:
            query_str = (f"xref_resolved:{b_strat}={b_xref}->CID{b_cid}"
                         f"({b_rt})")
            resolution_rows.append((
                cpd_id, "xref_resolved", query_str, "found",
                b_cid, b_ik, b_smi, time.time()))
            resolved_cpds.add(cpd_id)
            resolved_count += 1
        else:
            still_conflicting += 1
            still_conflicting_cpds.add(cpd_id)

    if resolution_rows:
        save_batch_to_cache(conn, db_lock, resolution_rows)

    logger.info("  Resolved: %d", resolved_count)
    logger.info("  Still conflicting: %d", still_conflicting)
    return resolved_cpds, still_conflicting_cpds


# ---------------------------------------------------------------------------
# Phase 2: Batched name lookups with multi-name cross-validation
# ---------------------------------------------------------------------------

def run_phase2_names(needs_name_lookup, names, structures, conn, db_lock,
                     workers):
    """Query PubChem by compound names with disambiguation.

    Returns (found_names, not_found_names, ambiguous_names) counts.
    """
    name_to_cpds = defaultdict(list)
    cpd_to_names = {}
    no_name = []

    for cpd_id in needs_name_lookup:
        cpd_names = names.get(cpd_id, [])
        best = pick_best_names(cpd_names)
        if best:
            cpd_to_names[cpd_id] = best
            for name in best:
                name_to_cpds[name].append(cpd_id)
        else:
            no_name.append(cpd_id)

    name_list = list(name_to_cpds.keys())
    logger.info("Phase 2: name lookups")
    logger.info("  Unique names to query: %d", len(name_list))
    logger.info("  Compounds with no usable name: %d", len(no_name))

    name_results = {}

    def do_name_batch(batch):
        return _query_names_batch_recursive(batch)

    batches = [name_list[i:i + NAME_BATCH_SIZE]
               for i in range(0, len(name_list), NAME_BATCH_SIZE)]

    if tqdm:
        pbar2 = tqdm(total=len(batches), desc="Phase 2: name batches")

    batches_done = 0
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(do_name_batch, batch): batch
                   for batch in batches}
        for future in as_completed(futures):
            batch_results = future.result()
            name_results.update(batch_results)
            batches_done += 1
            if batches_done % 10 == 0:
                logger.info("  Phase 2 progress: %d/%d batches done",
                            batches_done, len(batches))
            if tqdm:
                pbar2.update(1)

    if tqdm:
        pbar2.close()

    found_names = 0
    not_found_names = 0
    ambiguous_names = 0
    name_cache_rows = []

    for cpd_id, tried_names in cpd_to_names.items():
        found_cids = {}
        for name in tried_names:
            if name in name_results:
                cid, inchikey, smiles = name_results[name]
                found_cids[cid] = (name, inchikey, smiles)

        if not found_cids:
            name_cache_rows.append((cpd_id, "name_lookup", tried_names[0],
                                    "not_found", None, None, None, time.time()))
            not_found_names += 1
        elif len(found_cids) == 1:
            cid = list(found_cids.keys())[0]
            name, inchikey, smiles = found_cids[cid]
            name_cache_rows.append((cpd_id, "name_lookup", name, "found",
                                    cid, inchikey, smiles, time.time()))
            found_names += 1
        else:
            stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
            resolved = False
            if stored_ik:
                s_parts = stored_ik.split("-")
                # Tier 1: exact InChIKey match
                for cid, (name, inchikey, smiles) in found_cids.items():
                    if inchikey and inchikey == stored_ik:
                        name_cache_rows.append(
                            (cpd_id, "name_lookup", name, "found",
                             cid, inchikey, smiles, time.time()))
                        found_names += 1
                        resolved = True
                        break
                # Tier 2: blocks 1+2 match (connectivity + stereo)
                if not resolved:
                    for cid, (name, inchikey, smiles) in found_cids.items():
                        if not inchikey:
                            continue
                        p_parts = inchikey.split("-")
                        if (len(s_parts) >= 2 and len(p_parts) >= 2
                                and s_parts[0] == p_parts[0]
                                and s_parts[1] == p_parts[1]):
                            name_cache_rows.append(
                                (cpd_id, "name_lookup", name, "found",
                                 cid, inchikey, smiles, time.time()))
                            found_names += 1
                            resolved = True
                            break
                # Tier 3: block 1 only match (connectivity)
                if not resolved:
                    block1_matches = []
                    for cid, (name, inchikey, smiles) in found_cids.items():
                        if not inchikey:
                            continue
                        p_parts = inchikey.split("-")
                        if (len(s_parts) >= 1 and len(p_parts) >= 1
                                and s_parts[0] == p_parts[0]):
                            block1_matches.append(
                                (cid, name, inchikey, smiles))
                    if len(block1_matches) == 1:
                        cid, name, inchikey, smiles = block1_matches[0]
                        name_cache_rows.append(
                            (cpd_id, "name_lookup", name, "found",
                             cid, inchikey, smiles, time.time()))
                        found_names += 1
                        resolved = True
                    elif len(block1_matches) > 1:
                        stored_smi = structures.get(cpd_id, {}).get(
                            "smiles", "")
                        s_mol = (Chem.MolFromSmiles(stored_smi)
                                 if stored_smi else None)
                        if s_mol:
                            gen = rdFingerprintGenerator.GetMorganGenerator(
                                radius=2)
                            s_fp = gen.GetFingerprint(s_mol)
                            best_sim, best_match = -1, None
                            for cid, name, inchikey, smiles in block1_matches:
                                p_mol = (Chem.MolFromSmiles(smiles)
                                         if smiles else None)
                                if p_mol:
                                    p_fp = gen.GetFingerprint(p_mol)
                                    sim = DataStructs.TanimotoSimilarity(
                                        s_fp, p_fp)
                                    if sim > best_sim:
                                        best_sim = sim
                                        best_match = (cid, name, inchikey,
                                                      smiles)
                            if best_match:
                                cid, name, inchikey, smiles = best_match
                                name_cache_rows.append(
                                    (cpd_id, "name_lookup", name, "found",
                                     cid, inchikey, smiles, time.time()))
                                found_names += 1
                                resolved = True
            if not resolved:
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

    logger.info("  Found via name: %d", found_names)
    logger.info("  Not found: %d", not_found_names)
    if ambiguous_names:
        logger.info("  Ambiguous (names disagree): %d", ambiguous_names)

    return found_names, not_found_names, ambiguous_names


# ---------------------------------------------------------------------------
# Phase 2.5: InChIKey recovery for unresolved compounds
# ---------------------------------------------------------------------------

def run_phase25_inchikey_recovery(candidates, structures, conn, db_lock,
                                  workers):
    """Try direct InChIKey lookup for NOT_FOUND or MISMATCH compounds.

    Updates cache in-place for improved results.
    """
    cur_recovery = conn.execute(
        "SELECT cpd_id, strategy, status, pubchem_cid, pubchem_inchikey "
        "FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))),
        candidates,
    )
    recovery_candidates = []
    for row in cur_recovery.fetchall():
        cpd_id, strategy, status, pub_cid, pub_ik = row
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        if not stored_ik:
            continue
        if status == "not_found":
            recovery_candidates.append((cpd_id, stored_ik, "NOT_FOUND"))
        elif status == "found" and pub_ik:
            result_type = compare_inchikeys(stored_ik, pub_ik)
            if result_type in ("MISMATCH", "STEREO_DIFF"):
                recovery_candidates.append((cpd_id, stored_ik, result_type))

    ik_to_cpds = defaultdict(list)
    for cpd_id, stored_ik, old_rt in recovery_candidates:
        ik_to_cpds[stored_ik].append((cpd_id, old_rt))
    unique_iks = list(ik_to_cpds.keys())

    logger.info("Phase 2.5: InChIKey recovery")
    logger.info("  Recovery candidates: %d compounds (%d unique InChIKeys)",
                len(recovery_candidates), len(unique_iks))

    if unique_iks:
        ik_results = {}

        def do_ik_query(ik):
            return (ik, query_inchikey(ik))

        if tqdm:
            pbar_ik = tqdm(total=len(unique_iks),
                           desc="Phase 2.5: InChIKey lookups")

        ik_done = 0
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(do_ik_query, ik): ik
                       for ik in unique_iks}
            for future in as_completed(futures):
                ik, result = future.result()
                ik_results[ik] = result
                ik_done += 1
                if ik_done % 500 == 0:
                    logger.info("  Phase 2.5 progress: %d/%d InChIKey "
                                "lookups done", ik_done, len(unique_iks))
                if tqdm:
                    pbar_ik.update(1)

        if tqdm:
            pbar_ik.close()

        recovered = defaultdict(int)
        ik_cache_rows = []

        for stored_ik, cpd_list in ik_to_cpds.items():
            result = ik_results.get(stored_ik)
            if result is None:
                continue
            cid, pub_ik, smiles = result
            if not pub_ik:
                continue

            for cpd_id, old_rt in cpd_list:
                new_rt = compare_inchikeys(stored_ik, pub_ik)
                improved = False
                if old_rt == "NOT_FOUND":
                    improved = True
                elif old_rt == "MISMATCH" and new_rt in (
                        "MATCH", "PROTONATION_DIFF", "STEREO_DIFF"):
                    improved = True
                elif old_rt == "STEREO_DIFF" and new_rt in (
                        "MATCH", "PROTONATION_DIFF"):
                    improved = True

                if improved:
                    ik_cache_rows.append(
                        (cpd_id, "inchikey_lookup", stored_ik, "found",
                         cid, pub_ik, smiles, time.time()))
                    recovered[new_rt] += 1

        if ik_cache_rows:
            save_batch_to_cache(conn, db_lock, ik_cache_rows)

        total_recovered = sum(recovered.values())
        logger.info("  Recovered: %d", total_recovered)
        for rt in ["MATCH", "PROTONATION_DIFF", "STEREO_DIFF", "MISMATCH"]:
            if recovered[rt]:
                logger.info("    -> %s: %d", rt, recovered[rt])
        logger.info("  Still unresolved: %d",
                    len(recovery_candidates) - total_recovered)
    else:
        logger.info("  No candidates for InChIKey recovery.")


# ---------------------------------------------------------------------------
# Phase 2.6: InChI-based recovery for still-unresolved compounds
# ---------------------------------------------------------------------------

def run_phase26_inchi_recovery(candidates, structures, conn, db_lock,
                               workers):
    """Try direct InChI lookup for compounds still NOT_FOUND after Phase 2.5.

    Updates cache in-place for recovered results.
    """
    cur_inchi_recovery = conn.execute(
        "SELECT cpd_id, status FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))),
        candidates,
    )
    inchi_recovery_candidates = []
    for row in cur_inchi_recovery.fetchall():
        cpd_id, status = row
        if status == "not_found":
            stored_inchi = structures.get(cpd_id, {}).get("inchi", "")
            if stored_inchi:
                inchi_recovery_candidates.append((cpd_id, stored_inchi))

    logger.info("Phase 2.6: InChI-based recovery")
    logger.info("  Candidates with stored InChI: %d",
                len(inchi_recovery_candidates))

    if inchi_recovery_candidates:
        inchi_to_cpds = defaultdict(list)
        for cpd_id, stored_inchi in inchi_recovery_candidates:
            inchi_to_cpds[stored_inchi].append(cpd_id)
        unique_inchis = list(inchi_to_cpds.keys())

        def do_inchi_query(inchi_str):
            return (inchi_str, query_inchi(inchi_str))

        if tqdm:
            pbar_inchi = tqdm(total=len(unique_inchis),
                              desc="Phase 2.6: InChI lookups")

        inchi_results = {}
        inchi_done = 0
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(do_inchi_query, inchi): inchi
                       for inchi in unique_inchis}
            for future in as_completed(futures):
                inchi_str, result = future.result()
                inchi_results[inchi_str] = result
                inchi_done += 1
                if inchi_done % 500 == 0:
                    logger.info("  Phase 2.6 progress: %d/%d InChI "
                                "lookups done", inchi_done, len(unique_inchis))
                if tqdm:
                    pbar_inchi.update(1)

        if tqdm:
            pbar_inchi.close()

        inchi_recovered = 0
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
                inchi_recovered += 1

        if inchi_cache_rows:
            save_batch_to_cache(conn, db_lock, inchi_cache_rows)

        logger.info("  Recovered via InChI: %d", inchi_recovered)
    else:
        logger.info("  No candidates for InChI recovery.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Validate ModelSEED compound names against PubChem structures.")
    parser.add_argument("--limit", type=int, default=0,
                        help="Only process first N compounds (0 = all)")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from cache (default behavior)")
    parser.add_argument("--workers", type=int, default=10,
                        help="Number of concurrent request threads (default 10)")
    parser.add_argument("--apply", action="store_true",
                        help="Fetch full structures and apply corrections "
                             "for MISMATCH/STEREO_DIFF compounds")
    parser.add_argument("--output", type=str, default=None,
                        help="Output file for corrected structures "
                             "(default: writes to a new _corrected file, "
                             "never overwrites the original)")
    parser.add_argument("--images", action="store_true",
                        help="Generate side-by-side structure comparison "
                             "images in struct_imgs/")
    parser.add_argument("--max-images", type=int, default=20,
                        help="Max images per result type (default 20)")
    parser.add_argument("--clear-cache", action="store_true",
                        help="Clear the SQLite cache before running "
                             "(for a fresh re-validation)")
    parser.add_argument("--skip-ph7", action="store_true",
                        help="Skip pH 7 protonation normalization "
                             "when --apply is used")
    parser.add_argument("--review-report", action="store_true",
                        help="Generate detailed review report for "
                             "REVIEW:different_compound mismatches")
    args = parser.parse_args()

    # Handle cache clearing
    if args.clear_cache and os.path.exists(CACHE_DB):
        os.remove(CACHE_DB)
        logger.info("Cache cleared.")

    logger.info("Loading local data...")
    structures = load_structures(STRUCTURES_FILE)
    names = load_names(NAMES_FILE)
    external_ids = load_external_ids(ALIASES_FILE)

    # Build validation set: compounds with at least one name AND a stored InChIKey
    candidates = []
    for cpd_id in sorted(structures.keys()):
        if cpd_id in names and "inchikey" in structures[cpd_id]:
            candidates.append(cpd_id)

    logger.info("Total compounds with structures: %d", len(structures))
    logger.info("Total compounds with names: %d", len(names))
    logger.info("Validation candidates (name + InChIKey): %d", len(candidates))

    if args.limit > 0:
        candidates = candidates[:args.limit]
        logger.info("Limited to first %d compounds", args.limit)

    # Init cache
    conn = init_cache(CACHE_DB)
    db_lock = threading.Lock()
    cached = get_cached(conn)
    to_process = [c for c in candidates if c not in cached]
    logger.info("Already cached: %d", len(candidates) - len(to_process))
    logger.info("To process: %d", len(to_process))

    if not to_process:
        run_phase17_resolve_cached_conflicts(
            conn, db_lock, candidates, structures, workers=args.workers)
        generate_report(conn, candidates, structures, names=names,
                        report_file=REPORT_FILE, mismatch_file=MISMATCH_FILE,
                        protonation_file=PROTONATION_FILE,
                        stereo_file=STEREO_FILE)
        corrections = {}
        if args.apply or args.images:
            corrections = fetch_corrections(
                conn, db_lock, candidates, structures, workers=args.workers,
                names=names, rejected_file=REJECTED_FILE)
        # Normalize STEREO_DIFF corrections to pH 7 protonation
        if corrections and not args.skip_ph7:
            n = normalize_corrections_ph7(corrections)
            if n:
                logger.info("  STEREO_DIFF corrections normalized to pH 7: %d", n)
        # Phase 4: pH 7 protonation normalization
        if args.apply and not args.skip_ph7:
            ph7_corrections = run_phase4_protonation(
                candidates, structures, names=names,
                corrections_file=PROTONATION_CORRECTIONS_FILE)
            # Only add PH7 corrections for compounds not already corrected
            # (STEREO_DIFF corrections already have pH 7 normalization applied)
            for cpd_id, ph7_corr in ph7_corrections.items():
                if cpd_id not in corrections:
                    corrections[cpd_id] = ph7_corr
        if args.apply:
            output_path = args.output
            if output_path is None:
                output_path = os.path.join(
                    SCRIPT_DIR,
                    "Unique_ModelSEED_Structures_modified.txt")
            apply_corrections(corrections, structures,
                              output_path=output_path,
                              structures_file=STRUCTURES_FILE,
                              corrections_log=CORRECTIONS_LOG)
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
        conn.close()
        return

    # Phase 1: xref lookups
    (xref_resolved_cpds, xref_conflicts, xref_conflict_pairs,
     needs_name_lookup) = run_phase1_xref(
        to_process, external_ids, structures, conn, db_lock, args.workers)

    # Phase 1.5: resolve xref conflicts (in-memory data from Phase 1)
    resolved_cpds, _ = run_phase15_resolve_conflicts(
        xref_conflict_pairs, structures, conn, db_lock)
    xref_resolved_cpds.update(resolved_cpds)
    xref_conflicts -= resolved_cpds

    # Phase 1.7: resolve remaining xref conflicts from cache (fetches CID props)
    run_phase17_resolve_cached_conflicts(
        conn, db_lock, candidates, structures, workers=args.workers)

    # Compounds still needing name lookup
    chebi_to_cpds = defaultdict(list)
    kegg_to_cpds = defaultdict(list)
    for cpd_id in to_process:
        ext = external_ids.get(cpd_id, {"chebi": [], "kegg": []})
        for cid_val in ext["chebi"]:
            chebi_to_cpds[f"CHEBI:{cid_val}"].append(cpd_id)
        for kid in ext["kegg"]:
            kegg_to_cpds[kid].append(cpd_id)
    all_xref_cpds = set()
    for cpd_ids in chebi_to_cpds.values():
        all_xref_cpds.update(cpd_ids)
    for cpd_ids in kegg_to_cpds.values():
        all_xref_cpds.update(cpd_ids)
    unresolved_xref = [c for c in to_process if c in all_xref_cpds
                       and c not in xref_resolved_cpds
                       and c not in xref_conflicts]
    needs_name_lookup.extend(unresolved_xref)

    # Phase 2: name lookups
    run_phase2_names(needs_name_lookup, names, structures, conn, db_lock,
                     args.workers)

    # Phase 2.5: InChIKey recovery
    run_phase25_inchikey_recovery(candidates, structures, conn, db_lock,
                                  args.workers)

    # Phase 2.6: InChI recovery
    run_phase26_inchi_recovery(candidates, structures, conn, db_lock,
                               args.workers)

    # Generate report
    generate_report(conn, candidates, structures, names=names,
                    report_file=REPORT_FILE, mismatch_file=MISMATCH_FILE,
                    protonation_file=PROTONATION_FILE,
                    stereo_file=STEREO_FILE)

    # Phase 3: Fetch + apply corrections (only with --apply)
    corrections = {}
    if args.apply or args.images:
        corrections = fetch_corrections(
            conn, db_lock, candidates, structures, workers=args.workers,
            names=names, rejected_file=REJECTED_FILE)

    # Normalize STEREO_DIFF corrections to pH 7 protonation
    if corrections and not args.skip_ph7:
        n = normalize_corrections_ph7(corrections)
        if n:
            logger.info("  STEREO_DIFF corrections normalized to pH 7: %d", n)

    # Phase 4: pH 7 protonation normalization
    if args.apply and not args.skip_ph7:
        ph7_corrections = run_phase4_protonation(
            candidates, structures, names=names,
            corrections_file=PROTONATION_CORRECTIONS_FILE)
        # Only add PH7 corrections for compounds not already corrected
        # (STEREO_DIFF corrections already have pH 7 normalization applied)
        for cpd_id, ph7_corr in ph7_corrections.items():
            if cpd_id not in corrections:
                corrections[cpd_id] = ph7_corr

    if args.apply:
        output_path = args.output
        if output_path is None:
            output_path = os.path.join(
                SCRIPT_DIR,
                "Unique_ModelSEED_Structures_modified_v3.txt")
        apply_corrections(corrections, structures, output_path=output_path,
                          structures_file=STRUCTURES_FILE,
                          corrections_log=CORRECTIONS_LOG)

    if args.images:
        generate_comparison_images(conn, candidates, structures,
                                   images_dir=IMAGES_DIR,
                                   max_per_type=args.max_images,
                                   accepted_cpds=set(corrections.keys()))

    if args.review_report:
        generate_different_compound_review(
            conn, candidates, structures, names, external_ids,
            mismatch_file=MISMATCH_FILE, review_file=REVIEW_FILE,
            workers=args.workers)

    conn.close()


if __name__ == "__main__":
    main()
