#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

# Backfill the per-tool `pkas` dict onto every compound, recording each tool's
# pKa/pKb side by side rather than collapsing them into the single served
# pka/pkb fields, e.g.
#
#   "pkas": {
#       "Marvin":  {"pKa": "1:14:12.60;1:22:3.33;...",
#                   "pKb": "1:9:-3.03;1:14:-3.85;..."},
#       "MolGpKa": {"pKa": "1:14:10.12;1:29:0.98;...",
#                   "pKb": "1:6:3.66;1:9:3.42;..."}
#   }
#
# Each tool entry preserves the existing <microstate>:<atom>:<value> encoding
# and is keyed by the predicting tool (Marvin = the cascade-winning per-database
# Marvin value; MolGpKa = the OPAM2/MolGpKa prediction). A "" value means that
# tool predicted no ionizable atoms of that kind for the compound.
#
# The served top-level pka / pkb fields are NEVER touched: these per-tool records
# are added next to, not in place of, the picked values. The script only reads
# values already staged under Biochemistry/Structures/, so it can be re-run at
# any time and is idempotent.
#
# The Marvin selection mirrors Update_Compound_pKas.py exactly: only compounds
# with an accepted unique ModelSEED structure are considered, and the first DB
# in the priority cascade KEGG > MetaCyc > ChEBI > Rhea with a hit on one of the
# compound's aliases wins.

PKA_DBS = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]
MARVIN_LABEL = "Marvin"
MOLGPKA_LABEL = "MolGpKa"

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()

# --- Marvin per-database pKas, indexed by external id (cascade-winning value) ---
per_source_pkas = compounds_helper.loadPerSourcePkas(PKA_DBS)
cpd_pKab_dict = dict()
for (db, ext_id), entry in per_source_pkas.items():
    cpd_pKab_dict.setdefault(ext_id, dict())
    for kind, value in entry.items():
        cpd_pKab_dict[ext_id][kind] = value

# --- MolGpKa (OPAM2) pKas, keyed directly by ModelSEED compound id ---
OPAM2_PKA_FILE = os.path.dirname(__file__)+"/../../Biochemistry/Structures/ModelSEED/pkas/opam2_molgpka.tsv"
molgpka_pkas = dict()  # cpd -> {'pKa':str, 'pKb':str}
if(os.path.exists(OPAM2_PKA_FILE)):
    with open(OPAM2_PKA_FILE) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if(len(cols) < 3):
                continue
            molgpka_pkas.setdefault(cols[0], dict())[cols[1]] = cols[2]


def tool_entry(values):
    # Normalize a {'pKa':..,'pKb':..} mapping into an entry with both keys
    # present (empty string when absent). Returns None if neither is present.
    pka = values.get('pKa', "") if isinstance(values, dict) else ""
    pkb = values.get('pKb', "") if isinstance(values, dict) else ""
    if(pka in (None, "null")):
        pka = ""
    if(pkb in (None, "null")):
        pkb = ""
    if(pka == "" and pkb == ""):
        return None
    return {"pKa": pka, "pKb": pkb}


# Clear any pre-existing per-tool dict so a re-run cannot leave stale entries.
for cpd in compounds_dict:
    compounds_dict[cpd].pop('pkas', None)

marvin_count = 0
molgpka_count = 0
for cpd in sorted(compounds_dict.keys()):
    pkas = dict()

    # Marvin: only for compounds with an accepted unique ModelSEED structure,
    # taking the first DB in the cascade with a hit on one of its aliases.
    if(cpd in structures_dict):
        for DB in PKA_DBS:
            if(DB not in aliases_dict.get(cpd, dict())):
                continue
            hit = None
            for alias in aliases_dict[cpd][DB]:
                if(alias in cpd_pKab_dict):
                    hit = cpd_pKab_dict[alias]
                    break
            if(hit is not None):
                entry = tool_entry(hit)
                if(entry is not None):
                    pkas[MARVIN_LABEL] = entry
                    marvin_count += 1
                break

    # MolGpKa: direct compound-id lookup (also gated on a unique structure to
    # match how the served value is applied).
    if(cpd in structures_dict and cpd in molgpka_pkas):
        entry = tool_entry(molgpka_pkas[cpd])
        if(entry is not None):
            pkas[MOLGPKA_LABEL] = entry
            molgpka_count += 1

    # Additive: attach the per-tool dict next to the served pka/pkb, only for
    # compounds that actually have per-tool pKa data (compounds with none are
    # left without the field, matching how pka/pkb are empty when absent).
    if(pkas):
        compounds_dict[cpd]['pkas'] = pkas

print("Compounds with a Marvin pKa record:  "+str(marvin_count))
print("Compounds with a MolGpKa pKa record: "+str(molgpka_count))
print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
