import csv
import logging
import sqlite3
from collections import defaultdict

logger = logging.getLogger("pubchem_validate")


def load_structures(path):
    """Parse All_ModelSEED_Structures.txt -> {cpd_id: {smiles, inchikey, inchi, formula, charge}}.

    Format: cpd_id, type, charge_type, alias_id, source, formula, charge, structure
    (no header, 8 columns). Uses 'Charged' rows as canonical structures.
    """
    structs = {}
    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 8:
                continue
            cpd_id, typ, charge_type = row[0], row[1], row[2]
            formula, charge, structure = row[5], row[6], row[7]
            # Only use Charged rows as canonical
            if charge_type != "Charged":
                continue
            if cpd_id not in structs:
                structs[cpd_id] = {}
            # Store formula/charge from a row with non-null values
            if formula != "null" and "formula" not in structs[cpd_id]:
                structs[cpd_id]["formula"] = formula
                structs[cpd_id]["charge"] = charge
            key = typ.strip().lower()
            if structure and structure != "null":
                if key == "smile" and "smiles" not in structs[cpd_id]:
                    structs[cpd_id]["smiles"] = structure
                elif key == "inchikey" and "inchikey" not in structs[cpd_id]:
                    structs[cpd_id]["inchikey"] = structure
                elif key == "inchi" and "inchi" not in structs[cpd_id]:
                    structs[cpd_id]["inchi"] = structure
    return structs


def load_names(path):
    """Parse names file -> {cpd_id: [names]}."""
    names = defaultdict(list)
    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            if len(row) < 2:
                continue
            cpd_id, name = row[0], row[1]
            if name and name.strip():
                names[cpd_id].append(name.strip())
    return dict(names)


def load_external_ids(path):
    """Parse aliases file -> {cpd_id: {chebi: [...], kegg: [...]}}."""
    ids = defaultdict(lambda: {"chebi": [], "kegg": []})
    with open(path, "r") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)  # skip header
        for row in reader:
            if len(row) < 3:
                continue
            cpd_id, ext_id, source = row[0], row[1], row[2]
            if source == "ChEBI":
                ids[cpd_id]["chebi"].append(ext_id)
            elif source == "KEGG":
                ids[cpd_id]["kegg"].append(ext_id)
    return dict(ids)


def init_cache(db_path):
    """Create/open cache database, return connection."""
    conn = sqlite3.connect(db_path, check_same_thread=False)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS cache (
            cpd_id TEXT PRIMARY KEY,
            strategy TEXT,
            query TEXT,
            status TEXT,
            pubchem_cid TEXT,
            pubchem_inchikey TEXT,
            pubchem_smiles TEXT,
            timestamp REAL
        )
    """)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS corrections (
            cpd_id TEXT PRIMARY KEY,
            pubchem_cid TEXT,
            pubchem_smiles TEXT,
            pubchem_inchi TEXT,
            pubchem_inchikey TEXT,
            timestamp REAL
        )
    """)
    conn.commit()
    return conn


def get_cached(conn):
    """Return set of cpd_ids already in cache."""
    cur = conn.execute("SELECT cpd_id FROM cache")
    return {row[0] for row in cur.fetchall()}


def save_batch_to_cache(conn, lock, rows):
    """Write multiple rows to cache atomically."""
    with lock:
        conn.executemany(
            "INSERT OR REPLACE INTO cache VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.commit()
