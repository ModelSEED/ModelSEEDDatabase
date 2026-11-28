# deltag/deltagerr were removed from the records on 2026-09-11; energies now
# live per source under `thermodynamics`. Same precedence as the recommended
# direction. Sources with no energy are skipped: the LLMs entry carries a
# direction call and no number.
def _pick_energy(rec, missing=10000000):
    thermo = rec.get('thermodynamics')
    if isinstance(thermo, dict):
        for src in ('eQuilibrator', 'dGPredictor', 'Group contribution'):
            v = thermo.get(src)
            if isinstance(v, list) and len(v) >= 2:
                try:
                    return float(v[0]), float(v[1])
                except (TypeError, ValueError):
                    continue
    return float(missing), float(missing)

#!/usr/bin/env python
import os, sys
from BiochemPy import Compounds, Reactions,InChIs

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

print("\n================")
print("For Section: \"Computation of thermodynamic properties of ModelSEED compounds and reaction\"\n")

MS_Complete_Structures=dict()
with open("../../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt") as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')

        if("InChI" in array[5]):
            MS_Complete_Structures[array[5]]=1

# The MetaNetX chem_prop.tsv cross-check was retired (Biochemistry/Structures/
# MetaNetX/README.md): MetaNetX identifier versioning made it untenable, so
# the cache cross-check now lives in Structures_in_ModelSEED_and_eQuilibrator.txt.
# The InChI-overlap block is therefore omitted here.
print("InChI in ModelSEED: ",len(MS_Complete_Structures.keys()))
print("(MetaNetX cross-check omitted -- retired; see Biochemistry/Structures/MetaNetX/README.md)\n")

# The committed *_GroupContribution_eQuilibrator_Comparison.txt files are stale
# 2020-era artifacts no current script regenerates; compute the same stats from
# the live per-source `thermodynamics` dicts instead.
SOURCES = ("Group contribution", "eQuilibrator", "dGPredictor")
def _src(rec, label):
    v = (rec.get("thermodynamics") or {}).get(label)
    if isinstance(v, list) and len(v) >= 2:
        try:
            return float(v[0]), float(v[1])
        except (TypeError, ValueError):
            return None
    return None

Reported_Cpds = {"GC": 0, "EQ": 0, "DGP": 0, "EQR": 0, "LoEQE": 0, "HiEQE": 0, "GC_EQ": 0, "GC_DGP": 0, "EQ_DGP": 0}
for cpd, c in compounds_dict.items():
    if c.get("is_obsolete") == 1:
        continue
    gc, eq, dgp = _src(c, "Group contribution"), _src(c, "eQuilibrator"), _src(c, "dGPredictor")
    if gc and abs(gc[0]) < 1e6:
        Reported_Cpds["GC"] += 1
    if eq:
        Reported_Cpds["EQ"] += 1
        if eq[1] > 50: Reported_Cpds["EQR"] += 1
        if eq[1] < 5: Reported_Cpds["LoEQE"] += 1
        if eq[1] > 100: Reported_Cpds["HiEQE"] += 1
    if dgp and abs(dgp[0]) < 1e6:
        Reported_Cpds["DGP"] += 1
    if gc and eq: Reported_Cpds["GC_EQ"] += 1
    if gc and dgp: Reported_Cpds["GC_DGP"] += 1
    if eq and dgp: Reported_Cpds["EQ_DGP"] += 1

print("Compounds: ")
print("\tGC: ",Reported_Cpds["GC"])
print("\tEQ: ",Reported_Cpds["EQ"])
print("\tDGP: ",Reported_Cpds["DGP"])
print("\tRejected EQ: ",Reported_Cpds["EQR"])
print("\tShared GC&EQ: ",Reported_Cpds["GC_EQ"], " GC&DGP: ",Reported_Cpds["GC_DGP"], " EQ&DGP: ",Reported_Cpds["EQ_DGP"])

Reported_Rxns = {"GC": 0, "EQ": 0, "DGP": 0, "EQR": 0, "GCEQ": 0, "LoEQE": 0, "HiEQE": 0, "GC_DGP": 0, "EQ_DGP": 0}
for rxn, r in reactions_dict.items():
    if rxn not in reactions_dict or r.get("is_obsolete") == 1:
        continue
    gc, eq, dgp = _src(r, "Group contribution"), _src(r, "eQuilibrator"), _src(r, "dGPredictor")
    if gc and abs(gc[0]) < 1e6:
        Reported_Rxns["GC"] += 1
    if eq:
        Reported_Rxns["EQ"] += 1
        if gc and abs(gc[0]) < 1e6:
            Reported_Rxns["GCEQ"] += 1
        if eq[1] > 100: Reported_Rxns["EQR"] += 1
        if eq[1] < 5: Reported_Rxns["LoEQE"] += 1
        if eq[1] > 100: Reported_Rxns["HiEQE"] += 1
    if dgp and abs(dgp[0]) < 1e6:
        Reported_Rxns["DGP"] += 1
    if gc and dgp: Reported_Rxns["GC_DGP"] += 1
    if eq and dgp: Reported_Rxns["EQ_DGP"] += 1

print("Reactions: ")
print("\tGC: ",Reported_Rxns["GC"])
print("\tEQ: ",Reported_Rxns["EQ"])
print("\tDGP: ",Reported_Rxns["DGP"])
print("\tShared GC & EQ: ",Reported_Rxns["GCEQ"], " GC&DGP: ",Reported_Rxns["GC_DGP"], " EQ&DGP: ",Reported_Rxns["EQ_DGP"])
print("\tRejected EQ: ",Reported_Rxns["EQR"],"\n")

print("================")
print("For Section: \"Thermodynamics\"\n")

print("Compounds:")
pct = "{0:.2f}".format(float(Reported_Cpds["LoEQE"])/float(Reported_Cpds["EQ"]))
print("\tLow EQ error: ",Reported_Cpds["LoEQE"],pct)
pct = "{0:.2f}".format(float(Reported_Cpds["HiEQE"])/float(Reported_Cpds["EQ"]))
print("\tHigh EQ error: ",Reported_Cpds["HiEQE"],pct,"\n")
print("Reactions:")
pct = "{0:.2f}".format(float(Reported_Rxns["LoEQE"])/float(Reported_Rxns["EQ"]))
print("\tLow EQ error: ",Reported_Rxns["LoEQE"],pct)
pct = "{0:.2f}".format(float(Reported_Rxns["HiEQE"])/float(Reported_Rxns["EQ"]))
print("\tHigh EQ error: ",Reported_Rxns["HiEQE"],pct,"\n")

print("================")
print("For Table 4\n")

compound_counts=dict()
Cpd_Types = ['All', 'Total (GC)', 'Accepted (GC)', 'Total (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Cpd_Types:
    compound_counts[type]=0

for cpd in compounds_dict:
    if(compounds_dict[cpd]["is_obsolete"] == 1):
        continue

    cpd_obj = compounds_dict[cpd]

    compound_counts['All']+=1

    if('GC' in cpd_obj['notes']):
        compound_counts['Total (GC)']+=1

    if('EQ' in cpd_obj['notes']):
        compound_counts['Total (EQ)']+=1

    if('GC' in cpd_obj['notes'] or 'EQ' in cpd_obj['notes']):
        compound_counts['Structured']+=1

    if(_pick_energy(cpd_obj)[0] == 10000000):
        continue

    if('GC' in cpd_obj['notes'] and 'EQU' not in cpd_obj['notes']):
        compound_counts['Accepted (GC)']+=1
        compound_counts['Final']+=1
    elif('EQU' in cpd_obj['notes']):
        compound_counts['Accepted (EQ)']+=1
        compound_counts['Final']+=1

total_cpds = compound_counts['Structured']
print("Compounds Completeness")
for key in Cpd_Types:
    pct = "{0:.2f}".format(float(compound_counts[key])/float(total_cpds))
    print(key,compound_counts[key],pct)

reaction_counts=dict()
Rxn_Types = ['All', 'Total (GC)', 'Complete (GC)', 'Accepted (GC)', 'Total (EQ)', 'Complete (EQ)', 'Accepted (EQ)', 'Structured', 'Final']
for type in Rxn_Types:
    reaction_counts[type]=0

print("\nReactions Completeness")
for rxn in reactions_dict:
    if(reactions_dict[rxn]['status'] == 'EMPTY'):
        continue

    if(reactions_dict[rxn]["is_obsolete"] == 1):
        continue

    rxn_obj = reactions_dict[rxn]

    reaction_counts['All']+=1

    if('GCP' in rxn_obj['notes'] or 'GCC' in rxn_obj['notes']):
        reaction_counts['Total (GC)']+=1

    if('EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Total (EQ)']+=1

    if('GCP' in rxn_obj['notes'] or 'GCC' in rxn_obj['notes'] or \
           'EQP' in rxn_obj['notes'] or 'EQC' in rxn_obj['notes'] or 'Accepted (EQ)' in rxn_obj['notes']):
        reaction_counts['Structured']+=1

    if('GCC' in rxn_obj['notes']):
        reaction_counts['Complete (GC)']+=1

    if('EQC' in rxn_obj['notes']):
        reaction_counts['Complete (EQ)']+=1

    if(_pick_energy(rxn_obj)[0] == 10000000):
        continue

    if('EQU' in rxn_obj['notes']):
        reaction_counts['Accepted (EQ)']+=1
        reaction_counts['Final']+=1

    if('GCC' in rxn_obj['notes'] and 'EQU' not in rxn_obj['notes']):
        reaction_counts['Accepted (GC)']+=1
        reaction_counts['Final']+=1

total_rxns=reaction_counts['Structured']
for key in Rxn_Types:
    pct = "{0:.2f}".format(float(reaction_counts[key])/float(total_rxns))
    print(key,reaction_counts[key],pct)
print("\n================\n")
