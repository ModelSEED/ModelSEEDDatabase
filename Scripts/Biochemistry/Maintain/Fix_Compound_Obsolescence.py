#!/usr/bin/env python
from BiochemPy import Compounds

structs_check=['FKNQFGJONOIPTF-UHFFFAOYSA-N','AQLMHYSWFMLWBS-UHFFFAOYSA-N']

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["InChIKey"],["ModelSEED"])

InChIKey_Cpds=dict()
for cpd in sorted (structures_dict):
    if('InChIKey' not in structures_dict[cpd]):
        continue

    for structure in structures_dict[cpd]['InChIKey']:
        if(structure not in InChIKey_Cpds):
            InChIKey_Cpds[structure]=list()
        InChIKey_Cpds[structure].append(cpd)

Update_Compounds=0
for structure in sorted(InChIKey_Cpds):

    if(structure not in structs_check):
        continue

    primary_cpd = sorted(InChIKey_Cpds[structure])[0]

    if(len(InChIKey_Cpds[structure])==1):
        if(compounds_dict[primary_cpd]["is_obsolete"]==1):
            Update_Reactions+=1
            compounds_dict[primary_cpd]["is_obsolete"]=0
            compounds_dict[primary_cpd]["linked_reaction"]="null"
    else:
        for rxn in InChIKey_Cpds[structure]:
            cpd_list = ";".join(sorted(x for x in InChIKey_Cpds[structure] if x != rxn))
            if(rxn == primary_cpd and ( compounds_dict[rxn]["is_obsolete"]==1 or cpd_list != compounds_dict[rxn]["linked_compound"] )):
                print("Updating primary compound "+rxn+" and removing any indication that its obsolete")
                Update_Compounds+=1
                compounds_dict[rxn]["linked_compound"] = cpd_list
                compounds_dict[rxn]["is_obsolete"]=0
            elif(rxn != primary_cpd and ( compounds_dict[rxn]["is_obsolete"]==0 or cpd_list != compounds_dict[rxn]["linked_compound"] )):
                print("Updating compound "+rxn+" to indicate that its obsolete")
                Update_Compounds+=1
                compounds_dict[rxn]["linked_compound"] = cpd_list
                compounds_dict[rxn]["is_obsolete"]=1

if(Update_Compounds>0):
    print("Saving obsolescence updating in "+str(Update_Compounds)+" compounds")
    compounds_helper.saveCompounds(compounds_dict)



