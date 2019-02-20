#!/usr/bin/env python
import os, sys
temp=list();
header=1;

from BiochemPy import Compounds

CompoundsHelper = Compounds()
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
Compounds_Dict = CompoundsHelper.loadCompounds()

Structures_Root="../../Biochemistry/Structures/"
for cpd in sorted (Compounds_Dict.keys()):
    if(cpd not in Structures_Dict):
        Compounds_Dict[cpd]['inchikey']=""
        Compounds_Dict[cpd]['smiles']=""
    else:
        formula_charge_dict={'formula':'null','charge':'0'}
        inchikey=""
        smile=""
        if('InChI' in Structures_Dict[cpd]):
            struct = list(Structures_Dict[cpd]['InChI'].keys())[0]
            formula_charge_dict = Structures_Dict[cpd]['InChI'][struct]
        elif('SMILE' in Structures_Dict[cpd]):
            struct = list(Structures_Dict[cpd]['SMILE'].keys())[0]
            formula_charge_dict = Structures_Dict[cpd]['SMILE'][struct]
        if(formula_charge_dict['charge']=='null'):
            formula_charge_dict['charge']='0'

        if('InChIKey' in Structures_Dict[cpd]):
            inchikey = sorted(list(Structures_Dict[cpd]['InChIKey'].keys()))[0]
        if('SMILE' in Structures_Dict[cpd]):
            smile = sorted(list(Structures_Dict[cpd]['SMILE'].keys()))[0]

        Compounds_Dict[cpd]['inchikey']=inchikey
        Compounds_Dict[cpd]['smiles']=smile

        if(formula_charge_dict['formula'] != "null"):
            Compounds_Dict[cpd]['formula']=formula_charge_dict['formula']
            Compounds_Dict[cpd]['charge']=formula_charge_dict['charge']

print("Saving compounds")
CompoundsHelper.saveCompounds(Compounds_Dict)
