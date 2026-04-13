import csv
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import Pool, cpu_count

from rdkit import Chem

from pubchem_api import query_cid_properties
from structure_compare import (compare_inchikeys, compute_formula_charge_from_inchi,
                                validate_correction_worker)
from data_io import save_batch_to_cache

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, os.path.join(BASE_DIR, "ModelSEEDDatabase", "Libs", "Python"))

logger = logging.getLogger("pubchem_validate")


def fetch_corrections(conn, db_lock, candidates, structures, workers=5,
                      names=None, rejected_file=None):
    """Phase 3: fetch full structures for MISMATCH + STEREO_DIFF compounds.

    Applies validate_correction() gate to every candidate --- only returns
    corrections that demonstrably improve the data.

    Returns dict: {cpd_id: {pubchem_cid, smiles, inchi, inchikey, result_type,
                             strategy, query, validation_reason}}
    """
    # Identify correctable compounds from cache
    cur = conn.execute(
        "SELECT cpd_id, strategy, query, status, pubchem_cid, pubchem_inchikey "
        "FROM cache WHERE cpd_id IN ({})".format(
            ",".join("?" * len(candidates))),
        candidates,
    )
    correctable = {}  # cpd_id -> (cid, strategy, query, result_type)
    for row in cur.fetchall():
        cpd_id, strategy, query, status, pub_cid, pub_inchikey = row
        if status != "found" or not pub_cid:
            continue
        stored_ik = structures.get(cpd_id, {}).get("inchikey", "")
        result_type = compare_inchikeys(stored_ik, pub_inchikey)
        if result_type in ("MISMATCH", "STEREO_DIFF"):
            correctable[cpd_id] = (pub_cid, strategy, query, result_type)

    if not correctable:
        logger.info("Phase 3: No correctable compounds found.")
        return {}

    # Check which CIDs already have corrections cached
    cur2 = conn.execute("SELECT cpd_id FROM corrections")
    already_corrected = {row[0] for row in cur2.fetchall()}

    # Deduplicate by CID (multiple cpd_ids might map to same CID)
    unique_cids = {}  # cid -> [cpd_ids]
    for cpd_id, (cid, _, _, _) in correctable.items():
        if cpd_id not in already_corrected:
            unique_cids.setdefault(cid, []).append(cpd_id)

    logger.info("Phase 3: Fetching full structures for corrections")
    logger.info("  Correctable compounds: %d", len(correctable))
    logger.info("  Already cached corrections: %d",
                len(already_corrected & set(correctable)))
    logger.info("  Unique CIDs to fetch: %d", len(unique_cids))

    if unique_cids:
        if tqdm:
            pbar = tqdm(total=len(unique_cids), desc="Phase 3: CID lookups")

        def do_cid_query(cid):
            return (cid, query_cid_properties(cid))

        cid_done = 0
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(do_cid_query, cid): cid
                       for cid in unique_cids}
            for future in as_completed(futures):
                cid, props = future.result()
                cid_done += 1
                if cid_done % 100 == 0:
                    logger.info("  Phase 3 progress: %d/%d CID lookups done",
                                cid_done, len(unique_cids))
                if tqdm:
                    pbar.update(1)
                if props is None:
                    continue
                rows = []
                for cpd_id in unique_cids[cid]:
                    rows.append((cpd_id, cid, props["smiles"],
                                 props["inchi"], props["inchikey"], time.time()))
                with db_lock:
                    conn.executemany(
                        "INSERT OR REPLACE INTO corrections VALUES (?,?,?,?,?,?)",
                        rows,
                    )
                    conn.commit()

        if tqdm:
            pbar.close()

    # Build corrections dict from DB
    cur3 = conn.execute(
        "SELECT cpd_id, pubchem_cid, pubchem_smiles, pubchem_inchi, "
        "pubchem_inchikey FROM corrections WHERE cpd_id IN ({})".format(
            ",".join("?" * len(correctable))),
        list(correctable.keys()),
    )
    raw_result = {}
    for row in cur3.fetchall():
        cpd_id = row[0]
        cid, strategy, query, result_type = correctable[cpd_id]
        raw_result[cpd_id] = {
            "pubchem_cid": row[1],
            "smiles": row[2],
            "inchi": row[3],
            "inchikey": row[4],
            "result_type": result_type,
            "strategy": strategy,
            "query": query,
        }

    # Validate each correction: only accept those that improve the data
    # Use multiprocessing for the expensive validate_correction calls
    work_items = [(cpd_id, structures.get(cpd_id, {}), corr,
                   corr["result_type"])
                  for cpd_id, corr in raw_result.items()]

    validated = {}
    rejected = []
    num_workers = min(cpu_count(), 32)
    with Pool(num_workers) as pool:
        for cpd_id, accept, reason in pool.imap_unordered(
                validate_correction_worker, work_items, chunksize=32):
            corr = raw_result[cpd_id]
            if accept:
                corr["validation_reason"] = reason
                validated[cpd_id] = corr
            else:
                rejected.append((cpd_id, corr["result_type"], reason))

    logger.info("  Corrections validated: %d", len(validated))
    if rejected:
        logger.info("  Corrections rejected: %d", len(rejected))
        for cpd_id, rt, reason in rejected[:20]:
            logger.info("    %s (%s): %s", cpd_id, rt, reason)
        if len(rejected) > 20:
            logger.info("    ... and %d more", len(rejected) - 20)

        # Save rejected corrections to file for auditing
        if rejected_file is not None:
            with open(rejected_file, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t")
                writer.writerow(["cpd_id", "compound_name", "result_type", "reason",
                                 "pubchem_cid", "strategy", "query",
                                 "stored_smiles", "pubchem_smiles",
                                 "stored_inchikey", "pubchem_inchikey"])
                for cpd_id, rt, reason in rejected:
                    corr = raw_result.get(cpd_id, {})
                    stored = structures.get(cpd_id, {})
                    # Include compound name for easier review
                    cpd_names = names.get(cpd_id, []) if names else []
                    name_str = cpd_names[0] if cpd_names else ""
                    writer.writerow([
                        cpd_id, name_str, rt, reason,
                        corr.get("pubchem_cid", ""),
                        corr.get("strategy", ""),
                        corr.get("query", ""),
                        stored.get("smiles", ""),
                        corr.get("smiles", ""),
                        stored.get("inchikey", ""),
                        corr.get("inchikey", ""),
                    ])
            logger.info("  Rejected corrections: %s", rejected_file)

    return validated


def normalize_corrections_ph7(corrections):
    """Normalize STEREO_DIFF corrections to pH 7 protonation.

    When a STEREO_DIFF correction replaces the stored InChI with PubChem's
    version, the protonation state changes from the curated ionized form
    (appropriate for cellular pH ~7) to PubChem's neutral form.  This function
    fixes that by running OpenBabel CorrectForPH(7.0) on the PubChem SMILES
    and updating the correction in-place when only protonation changed.

    Returns the number of corrections that were pH 7-normalized.
    """
    from protonation import normalize_to_ph7

    stereo_cpds = [(cpd_id, corr) for cpd_id, corr in corrections.items()
                   if corr.get("result_type") == "STEREO_DIFF"]

    if not stereo_cpds:
        return 0

    normalized = 0
    for cpd_id, corr in stereo_cpds:
        pub_smiles = corr.get("smiles", "")
        if not pub_smiles:
            continue

        ph7 = normalize_to_ph7(pub_smiles)
        if ph7 is None:
            continue

        # Compare PubChem InChIKey to pH 7 InChIKey
        pub_ik = corr.get("inchikey", "")
        ph7_ik = ph7["inchikey"]
        result = compare_inchikeys(pub_ik, ph7_ik)

        if result == "MATCH":
            # Already at pH 7, no change needed
            continue

        if result == "PROTONATION_DIFF":
            # pH 7 changes only protonation — update correction in-place
            corr["smiles"] = ph7["smiles"]
            corr["inchi"] = ph7["inchi"]
            corr["inchikey"] = ph7["inchikey"]
            normalized += 1
            logger.info("  %s: STEREO_DIFF correction normalized to pH 7 "
                        "(InChIKey %s -> %s)", cpd_id,
                        pub_ik[-1], ph7_ik[-1])

        # If STEREO_DIFF or MISMATCH: OpenBabel changed connectivity/stereo
        # — keep PubChem's version as-is (normalization artifact)

    return normalized


def apply_corrections(corrections, structures, output_path=None,
                      structures_file=None, corrections_log=None):
    """Apply corrections and write to a NEW output file in Unique format.

    Reads All_ModelSEED_Structures.txt (8-col, no header), keeps only Charged
    rows, merges aliases per (cpd_id, type), applies corrections, and writes
    in Unique format: ID, Type, Aliases, Formula, Charge, Structure (with header).
    Never overwrites the original All file.
    """
    if structures_file is None:
        structures_file = os.path.join(SCRIPT_DIR, "All_ModelSEED_Structures.txt")
    if corrections_log is None:
        corrections_log = os.path.join(SCRIPT_DIR, "pubchem_corrections_log.tsv")

    # Determine output path (never overwrite original)
    if output_path is None:
        output_path = os.path.join(SCRIPT_DIR, "Unique_ModelSEED_Structures_modified.txt")
    logger.info("Original file: %s (unchanged)", structures_file)
    logger.info("Output file:   %s", output_path)

    # Read All file and build Unique representation (Charged only, merged aliases)
    # Key: (cpd_id, type) -> {aliases: set, formula, charge, structure}
    unique = {}  # (cpd_id, typ) -> dict
    cpd_order = []  # preserve compound ordering
    seen_cpds = set()
    with open(structures_file, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 8:
                continue
            cpd_id, typ, charge_type = row[0], row[1], row[2]
            alias_id, formula, charge, structure = row[3], row[5], row[6], row[7]
            if charge_type != "Charged":
                continue
            key = (cpd_id, typ)
            if key not in unique:
                unique[key] = {
                    "aliases": set(),
                    "formula": formula if formula != "null" else "",
                    "charge": charge if charge != "null" else "",
                    "structure": structure if structure != "null" else "",
                }
            unique[key]["aliases"].add(alias_id)
            # Update formula/charge if current is empty and this row has values
            if not unique[key]["formula"] and formula != "null":
                unique[key]["formula"] = formula
                unique[key]["charge"] = charge if charge != "null" else ""
            if cpd_id not in seen_cpds:
                cpd_order.append(cpd_id)
                seen_cpds.add(cpd_id)

    # Propagate formula/charge across types for the same compound
    # (InChIKey rows in the All file often have "null" formula/charge)
    for cpd_id in cpd_order:
        # Find the best formula/charge from any type for this compound
        best_formula, best_charge = "", ""
        for typ in ["SMILE", "InChI", "InChIKey"]:
            key = (cpd_id, typ)
            if key in unique and unique[key]["formula"]:
                best_formula = unique[key]["formula"]
                best_charge = unique[key]["charge"]
                break
        # Apply to all types that are missing formula/charge
        for typ in ["SMILE", "InChI", "InChIKey"]:
            key = (cpd_id, typ)
            if key in unique and not unique[key]["formula"]:
                unique[key]["formula"] = best_formula
                unique[key]["charge"] = best_charge

    # Compute new formula/charge for corrected compounds
    formula_updates = {}  # cpd_id -> (new_formula, new_charge)
    for cpd_id, corr in corrections.items():
        new_inchi = corr.get("inchi", "")
        if new_inchi:
            new_formula, new_charge = compute_formula_charge_from_inchi(new_inchi)
            if new_formula is not None:
                formula_updates[cpd_id] = (
                    new_formula,
                    str(new_charge) if new_charge is not None else "0"
                )

    # Open corrections log for incremental writing
    log_header = [
        "timestamp", "cpd_id", "result_type", "field",
        "old_value", "new_value",
        "pubchem_cid", "strategy", "query",
    ]
    # Overwrite (not append) so the log matches the output file, which is
    # also regenerated from scratch each run.
    log_fh = open(corrections_log, "w", newline="")
    log_writer = csv.writer(log_fh, delimiter="\t")
    log_writer.writerow(log_header)

    # Apply corrections to the unique representation, logging each change
    total_changes = 0
    corrected_cpds = set()
    type_to_field = {"SMILE": "smiles", "InChIKey": "inchikey", "InChI": "inchi"}
    ts = datetime.now().isoformat()

    for (cpd_id, typ), entry in unique.items():
        if cpd_id not in corrections:
            continue
        corr = corrections[cpd_id]

        # Update structure value
        field_key = type_to_field.get(typ)
        if field_key is not None:
            old_val = entry["structure"]
            new_val = corr[field_key]
            if new_val and old_val != new_val:
                entry["structure"] = new_val
                log_writer.writerow([ts, cpd_id, corr["result_type"], typ,
                                     old_val, new_val, corr["pubchem_cid"],
                                     corr["strategy"], corr["query"]])
                log_fh.flush()
                total_changes += 1
                corrected_cpds.add(cpd_id)

        # Update formula/charge
        if cpd_id in formula_updates:
            new_formula, new_charge = formula_updates[cpd_id]
            entry["formula"] = new_formula
            entry["charge"] = new_charge

    # Log formula/charge changes (once per compound)
    for cpd_id, (new_formula, new_charge) in formula_updates.items():
        old_formula = structures.get(cpd_id, {}).get("formula", "")
        old_charge = structures.get(cpd_id, {}).get("charge", "")
        corr = corrections[cpd_id]
        if old_formula != new_formula or old_charge != new_charge:
            log_writer.writerow([ts, cpd_id, corr["result_type"],
                                 "Formula/Charge",
                                 f"{old_formula}/{old_charge}",
                                 f"{new_formula}/{new_charge}",
                                 corr["pubchem_cid"], corr["strategy"],
                                 corr["query"]])
            log_fh.flush()
            total_changes += 1
            corrected_cpds.add(cpd_id)

    log_fh.close()

    if not corrections:
        logger.info("No corrections to apply.")

    # Write in Unique format: ID, Type, Aliases, Formula, Charge, Structure
    # Use atomic write (temp file + rename) for crash safety
    type_order = ["SMILE", "InChIKey", "InChI"]
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["ID", "Type", "Aliases", "Formula", "Charge", "Structure"])
        for cpd_id in cpd_order:
            for typ in type_order:
                key = (cpd_id, typ)
                if key not in unique:
                    continue
                entry = unique[key]
                aliases_str = ";".join(sorted(entry["aliases"]))
                writer.writerow([
                    cpd_id, typ, aliases_str,
                    entry["formula"], entry["charge"], entry["structure"],
                ])
    os.replace(tmp_path, output_path)

    logger.info("Output written in Unique format (Charged only, aliases merged).")
    logger.info("  Compounds corrected: %d", len(corrected_cpds))
    logger.info("  Fields changed: %d", total_changes)
    logger.info("  Output file:     %s", output_path)
    logger.info("  Corrections log: %s", corrections_log)

    # Post-correction self-consistency check on corrected compounds
    if corrected_cpds:
        logger.info("  Running post-correction consistency checks...")
        inconsistencies = 0
        for cpd_id in sorted(corrected_cpds):
            cpd_data = {}
            for typ in type_order:
                key = (cpd_id, typ)
                if key in unique:
                    cpd_data[typ] = unique[key]["structure"]
            smiles = cpd_data.get("SMILE", "")
            inchi = cpd_data.get("InChI", "")
            inchikey = cpd_data.get("InChIKey", "")
            issues = []
            # Check SMILES parses
            smi_mol = Chem.MolFromSmiles(smiles) if smiles else None
            if smiles and smi_mol is None:
                issues.append("SMILES fails RDKit parsing")
            # Check InChI parses
            inchi_mol = Chem.inchi.MolFromInchi(inchi) if inchi else None
            if inchi and inchi_mol is None:
                issues.append("InChI fails RDKit parsing")
            # Verify InChIKey matches InChI
            if inchi and inchikey:
                try:
                    computed_key = Chem.inchi.InchiToInchiKey(inchi)
                    if computed_key and computed_key != inchikey:
                        issues.append(
                            f"InChIKey mismatch: stored={inchikey[:14]} "
                            f"computed={computed_key[:14]}")
                except Exception:
                    issues.append("InChIKey computation failed")
            # Cross-check SMILES vs InChI
            if smi_mol and inchi_mol:
                try:
                    key_smi = Chem.inchi.MolToInchiKey(smi_mol)
                    key_inchi = Chem.inchi.InchiToInchiKey(inchi)
                    if key_smi and key_inchi:
                        if key_smi.split("-")[0] != key_inchi.split("-")[0]:
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
