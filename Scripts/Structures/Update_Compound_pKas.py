#!/usr/bin/env python
#
# Picks the single served pka/pkb per compound from the Marvin per-database
# cascade (KEGG > MetaCyc > ChEBI > Rhea). The companion script
# Add_Compound_pKa_Sources.py then records each tool's pKa/pKb ADDITIVELY in a
# per-compound `pkas` dict (Marvin = this cascade-winning value; MolGpKa = the
# OPAM2/MolGpKa prediction), next to these served values without changing them;
# the cascade runs it as the AddPkaSources stage right after this one.
import os,sys
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()

# Load pKas and pKbs from the post-A1 layout: <db>/pkas/<tool>_<ver>.tsv.
# Iteration order is the priority cascade KEGG > MetaCyc > ChEBI > Rhea.
# First DB with a hit on any of the compound's aliases wins.
#
# ChEBI and Rhea were previously unused because of ID-format mismatches
# (ChEBI: 'CHEBI_15377' in pKa file vs '15377' in aliases; Rhea:
# 'POLYMER_X' vs 'POLYMER:X'). Compounds.loadPerSourcePkas normalizes
# both at ingestion now, see sources.yaml for the migration history.
# Adding ChEBI to the cascade unlocks ~3,800 additional pKa
# attributions; Rhea's polymer pKa data is fully shadowed by primary
# sources so its inclusion is for completeness rather than coverage.
PKA_DBS = ["KEGG", "MetaCyc", "ChEBI", "Rhea"]
per_source_pkas = compounds_helper.loadPerSourcePkas(PKA_DBS)
cpd_pKab_dict = dict()
for (db, ext_id), entry in per_source_pkas.items():
    if(ext_id not in cpd_pKab_dict):
        cpd_pKab_dict[ext_id] = dict()
    for kind, value in entry.items():
        cpd_pKab_dict[ext_id][kind] = value

# We're removing all pKa and pKb before loading new ones
for cpd in compounds_dict:
    compounds_dict[cpd]['pka']=""
    compounds_dict[cpd]['pkb']=""

# We're only loading pKa/pKb for compounds that have an accepted unique structure in ModelSEED
for cpd in structures_dict:
    found=False
    for DB in PKA_DBS:
        if(found is True or DB not in aliases_dict[cpd]):
            continue

        for alias in aliases_dict[cpd][DB]:
            if(alias in cpd_pKab_dict):
                print(cpd,alias,cpd_pKab_dict[alias])
                if('pKa' in cpd_pKab_dict[alias]):
                    print(cpd,alias)
                    compounds_dict[cpd]['pka']=cpd_pKab_dict[alias]['pKa']
                else:
                    compounds_dict[cpd]['pka']=""

                if('pKb' in cpd_pKab_dict[alias]):
                    compounds_dict[cpd]['pkb']=cpd_pKab_dict[alias]['pKb']
                else:
                    compounds_dict[cpd]['pkb']=""

                # All structures for the same compound should be the same,
                # and so only need to process once
                found=True
                break

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
