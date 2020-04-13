#!/usr/bin/env python
import os
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
compounds_dict = compounds_helper.loadCompounds()

inchikey_dict = dict()
for cpd in structures_dict:
    if('InChIKey' in structures_dict[cpd]):
        for struct in structures_dict[cpd]['InChIKey'].keys():
            if(struct not in inchikey_dict):
                inchikey_dict[struct]=list()
            inchikey_dict[struct].append(cpd)

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
for cpd in sorted (compounds_dict.keys()):
    if(cpd not in structures_dict):
        compounds_dict[cpd]['inchikey']=""
        compounds_dict[cpd]['smiles']=""

        formula_charge_dict={'formula':'null','charge':0}
        #Override manually
        if(cpd in Overridden_Fields and 'formula' in Overridden_Fields[cpd]):
            formula_charge_dict['formula']=Overridden_Fields[cpd]['formula']
        if(cpd in Overridden_Fields and 'charge' in Overridden_Fields[cpd]):
            formula_charge_dict['charge']=int(Overridden_Fields[cpd]['charge'])

        if(formula_charge_dict['formula'] != "null"):
            compounds_dict[cpd]['formula']=formula_charge_dict['formula']
            compounds_dict[cpd]['charge']=int(formula_charge_dict['charge'])
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
