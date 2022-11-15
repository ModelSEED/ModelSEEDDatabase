#!/usr/bin/env python
import os,sys
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

Overridden_Fields=dict()
header=list()
with open(os.path.dirname(__file__)+'/ACPs_Master_Formula_Charge.txt') as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')

        cpd=array.pop(0)

        if(len(header)==0):
            header=array
            continue

        if(cpd not in Overridden_Fields):
            Overridden_Fields[cpd]=dict()

        for i in range(len(array)):
            if(array[i] == 'null' or array[i] == 10000000):
                continue
            Overridden_Fields[cpd][header[i]]=array[i]

compounds_helper = Compounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
ignoring_structures_dict = compounds_helper.loadStructures(["SMILE","InChIKey"],["KEGG","MetaCyc"])
aliases_dict = compounds_helper.loadSourceAliases()
compounds_dict = compounds_helper.loadCompounds()
Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"

inchikey_dict = dict()
smiles_dict = dict()
for cpd in structures_dict:
    if('InChIKey' in structures_dict[cpd]):
        for struct in structures_dict[cpd]['InChIKey'].keys():
            if(struct not in inchikey_dict):
                inchikey_dict[struct]=list()
            inchikey_dict[struct].append(cpd)
    if('SMILE' in structures_dict[cpd]):
        for struct in structures_dict[cpd]['SMILE'].keys():
            if(struct not in smiles_dict):
                smiles_dict[struct]=list()
            smiles_dict[struct].append(cpd)

Ignored_Structures=list()
with open(Structures_Root+"Curation/Ignored_Structures_Publication2020.txt") as fh:
    for line in fh.readlines():
        line=line.strip()
        array=line.split('\t')
        if(array[1] != "None"):
            continue

        external_id = array.pop(0)

        if(external_id in ignoring_structures_dict['SMILE']):
            if(len(aliases_dict['MetaCyc'][external_id])==1):
                cpd = list(aliases_dict['MetaCyc'][external_id])[0]
                if(cpd not in Ignored_Structures):
                    Ignored_Structures.append(cpd)
            else:
                print("Warning, multiple MS compounds for "+external_id)
        else:
            print("Warning, no SMILES to ignore for "+external_id)

for cpd in sorted (compounds_dict.keys()):
    if(cpd not in structures_dict):
        compounds_dict[cpd]['inchikey']=""
        compounds_dict[cpd]['smiles']=""

        if(cpd in Ignored_Structures):
            compounds_dict[cpd]['formula']="R"
            compounds_dict[cpd]['charge']="0"
            continue

        if(cpd in Overridden_Fields):
            if('formula' in Overridden_Fields[cpd]):
                compounds_dict[cpd]['formula']=Overridden_Fields[cpd]['formula']
            if('charge' in Overridden_Fields[cpd]):
                compounds_dict[cpd]['charge']=int(Overridden_Fields[cpd]['charge'])
            continue

    else:
        formula_charge_dict={'formula':'null','charge':'0'}
        inchikey=""
        smile=""
        if('InChI' in structures_dict[cpd]):
            struct = list(structures_dict[cpd]['InChI'].keys())[0]
            formula_charge_dict = structures_dict[cpd]['InChI'][struct]
        elif('SMILE' in structures_dict[cpd]):
            struct = list(structures_dict[cpd]['SMILE'].keys())[0]
            formula_charge_dict = structures_dict[cpd]['SMILE'][struct]
        if(formula_charge_dict['charge']=='null'):
            formula_charge_dict['charge']='0'

        if('InChIKey' in structures_dict[cpd]):
            inchikey = sorted(list(structures_dict[cpd]['InChIKey'].keys()))[0]
        if('SMILE' in structures_dict[cpd]):
            smile = sorted(list(structures_dict[cpd]['SMILE'].keys()))[0]

        #See if any duplicates for newly assigned InChIKeys
        if(compounds_dict[cpd]['inchikey']==""):
            if(inchikey in inchikey_dict and len(inchikey_dict[inchikey])>1):
                print("Warning: Duplicate InChIKey: "+inchikey+" in "+" and ".join(inchikey_dict[inchikey]))
        compounds_dict[cpd]['inchikey']=inchikey
        compounds_dict[cpd]['smiles']=smile

        if(formula_charge_dict['formula'] != "null"):
            compounds_dict[cpd]['formula']=formula_charge_dict['formula']
            compounds_dict[cpd]['charge']=int(formula_charge_dict['charge'])

        #Override manually
        if(cpd in Overridden_Fields and 'formula' in Overridden_Fields[cpd]):
            compounds_dict[cpd]['formula']=Overridden_Fields[cpd]['formula']
        if(cpd in Overridden_Fields and 'charge' in Overridden_Fields[cpd]):
            compounds_dict[cpd]['charge']=int(Overridden_Fields[cpd]['charge'])

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
