#!/usr/bin/env python

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import os,sys,re
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

# REWRITTEN 2026-09-07. The previous version had four defects, all of which
# only bit when it was re-run -- which had not happened since the MolGpKa
# refresh and the Literature addition landed:
#
#   1. it hardcoded pkas/molgpka_opam2.tsv, the superseded 22,395-compound
#      export, while sources.yaml flags molgpka_stock-2026.09.tsv (32,201) as
#      consumed_by_production. Re-running downgraded MolGpKa by 9,747 records.
#   2. it cleared `pkas` and rebuilt it from Marvin and MolGpKa only, so the 31
#      Alberty-derived Literature entries were deleted every run.
#   3. it omitted the `kind` key, which pka_encoding.ladder() requires to refuse
#      building a sequential ladder out of per-site values.
#   4. it copied the source tables' THREE-field <fragment>:<atom>:<value>
#      strings verbatim, which pka_encoding.decode() rejects by design. The
#      well-formed two-field data in the JSONs came from the A1 migration; this
#      script would have undone it.
#
# It is now manifest-driven: which MolGpKa table is production is a sources.yaml
# edit, not a code edit, matching how the eQuilibrator pipeline resolves it.

PKA_DBS = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]
MARVIN_LABEL, MOLGPKA_LABEL, LITERATURE_LABEL = "Marvin", "MolGpKa", "Literature"

HERE = os.path.dirname(os.path.abspath(__file__))
PKA_DIR = os.path.join(HERE, "..", "..", "Biochemistry", "Structures", "ModelSEED", "pkas")
MANIFEST = os.path.join(HERE, "..", "..", "Biochemistry", "Structures", "sources.yaml")

sys.path.append(HERE)
from pka_encoding import migrate, MICROSCOPIC, MACROSCOPIC


def production_molgpka_table():
    """The MolGpKa table sources.yaml flags consumed_by_production.

    Parsed without a YAML dependency -- this script runs in the plain cascade
    environment. Falls back to the stock table by name, and only then to the
    legacy opam2 export, so a manifest typo degrades loudly rather than
    silently reinstating 22,395-compound coverage.
    """
    want = None
    try:
        with open(MANIFEST) as fh:
            current = None
            for line in fh:
                m = re.search(r"file:\s+(pkas/molgpka_[^\s]+)", line)
                if m:
                    current = m.group(1)
                elif current and "consumed_by_production:" in line:
                    if line.split(":", 1)[1].strip().lower() == "true":
                        want = current
                    current = None
    except OSError:
        pass
    for cand in (want, "pkas/molgpka_stock-2026.09.tsv", "pkas/molgpka_opam2.tsv"):
        if not cand:
            continue
        path = os.path.join(PKA_DIR, os.path.basename(cand))
        if os.path.exists(path):
            return path
    return None


def read_seed_keyed(path, id_col=0, kind_col=1, val_col=2):
    """<seed_id, kind, value> table -> {cpd: {'pKa':.., 'pKb':..}}."""
    out = dict()
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(id_col, kind_col, val_col):
                continue
            out.setdefault(cols[id_col], dict())[cols[kind_col]] = cols[val_col]
    return out


compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE", "InChI", "InChIKey"], ["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()

per_source_pkas = compounds_helper.loadPerSourcePkas(PKA_DBS)
cpd_pKab_dict = dict()
for (db, ext_id), entry in per_source_pkas.items():
    cpd_pKab_dict.setdefault(ext_id, dict())
    for kind, value in entry.items():
        cpd_pKab_dict[ext_id][kind] = value

molgpka_table = production_molgpka_table()
molgpka_pkas = read_seed_keyed(molgpka_table)
print("MolGpKa table:   " + (os.path.basename(molgpka_table) if molgpka_table else "NONE FOUND"))

# LITERATURE IS PRESERVED, NOT REBUILT. Update_Compound_Literature_pKas.py owns
# that entry: it reads Alberty's BasicBiochemData3 package directly and yields
# 31 compounds. The redistributable table in pkas/ is a ranked subset (8
# compounds) written by build_ladder_table.py for the evidence work, so
# rebuilding Literature from it would silently drop 23. Carrying the existing
# entry forward keeps this script idempotent and keeps one writer per source.
existing_literature = {
    cpd: rec["pkas"][LITERATURE_LABEL]
    for cpd, rec in compounds_dict.items()
    if isinstance(rec.get("pkas"), dict) and LITERATURE_LABEL in rec["pkas"]
}
print("Literature entries carried forward: " + str(len(existing_literature)))


def tool_entry(values, kind):
    """Normalise to {'kind', 'pKa', 'pKb'} in the two-field encoding.

    migrate() drops the atom index; it is idempotent, so a table already in the
    two-field form passes through unchanged.
    """
    pka = values.get("pKa", "") if isinstance(values, dict) else ""
    pkb = values.get("pKb", "") if isinstance(values, dict) else ""
    pka = "" if pka in (None, "null") else migrate(pka)
    pkb = "" if pkb in (None, "null") else migrate(pkb)
    if pka == "" and pkb == "":
        return None
    return {"kind": kind, "pKa": pka, "pKb": pkb}


for cpd in compounds_dict:
    compounds_dict[cpd].pop("pkas", None)

counts = {MARVIN_LABEL: 0, MOLGPKA_LABEL: 0, LITERATURE_LABEL: 0}
for cpd in sorted(compounds_dict.keys()):
    pkas = dict()
    if cpd in structures_dict:
        for DB in PKA_DBS:
            if DB not in aliases_dict.get(cpd, dict()):
                continue
            hit = None
            for alias in aliases_dict[cpd][DB]:
                if alias in cpd_pKab_dict:
                    hit = cpd_pKab_dict[alias]
                    break
            if hit is not None:
                entry = tool_entry(hit, MICROSCOPIC)
                if entry is not None:
                    pkas[MARVIN_LABEL] = entry
                    counts[MARVIN_LABEL] += 1
                break
    if cpd in structures_dict and cpd in molgpka_pkas:
        entry = tool_entry(molgpka_pkas[cpd], MICROSCOPIC)
        if entry is not None:
            pkas[MOLGPKA_LABEL] = entry
            counts[MOLGPKA_LABEL] += 1
    if cpd in existing_literature:
        pkas[LITERATURE_LABEL] = existing_literature[cpd]
        counts[LITERATURE_LABEL] += 1
    if pkas:
        compounds_dict[cpd]["pkas"] = pkas

for label in (MARVIN_LABEL, MOLGPKA_LABEL, LITERATURE_LABEL):
    print("Compounds with a %-10s pKa record: %d" % (label, counts[label]))
print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
