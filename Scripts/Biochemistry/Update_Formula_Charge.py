#!/usr/bin/env python
import os, sys, re
temp=list();
header=1;

import pybel
from rdkit.Chem import AllChem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Structures_Dict = CompoundsHelper.loadStructures(["InChI","SMILE"],["ModelSEED"])

Update_Compounds=0
Update_Formula=0
Update_Charge=0
for cpd in sorted(Compounds_Dict.keys()):
    if(cpd not in Structures_Dict):
        continue

    print(cpd)

    if('InChI' not in Structures_Dict[cpd] and 'SMILE' not in Structures_Dict[cpd]):
        continue

    print(cpd)
    
    mol=None
    mol_source=""
    try:
        if('InChI' in Structures_Dict[cpd]):
            mol = AllChem.MolFromInchi(Structures_Dict[cpd]['InChI'])
            if(mol==None):
                mol=pybel.readstring("inchi",Structures_Dict[cpd]['InChI'])
                if(mol):
                    mol_source="OpenBabel"
            else:
                mol_source="RDKit"

        elif('SMILE' in Structures_Dict[cpd]):
            mol = AllChem.MolFromSmiles(Structures_Dict[cpd]['SMILE'])
            if(mol==None):
                mol=pybel.readstring("smiles",Structures_Dict[cpd]['SMILE'])
                if(mol):
                    mol_source="OpenBabel"
            else:
                mol_source="RDKit"
    except Exception as e:
        pass

    if(mol == None):
        continue

    new_formula=""
    if(mol_source=="RDKit"):
        new_formula = AllChem.CalcMolFormula(mol)
    elif(mol_source=="OpenBabel"):
        new_formula=mol.formula

    new_charge=0
    if(mol_source=="RDKit"):
        new_charge = AllChem.GetFormalCharge(mol)
        match = re.search('([-+]\d?)$',new_formula)
        if(match):
            new_formula = new_formula.replace(match.group(),'')
    elif(mol_source=="OpenBabel"):
        new_charge = mol.charge
        match = re.search('([-+]+)$',new_formula)
        if(match):
            new_formula = new_formula.replace(match.group(),'')

    #normalizing formula my own way, so I can be consistent
    #these are hill-sorted, and merges molecular fragments
    new_formula = Compounds.mergeFormula(new_formula)[0]

    current_formula = Compounds_Dict[cpd]['formula']
    if(new_formula != current_formula):
        print("Updating formula for "+cpd+" with "+new_formula)
        Compounds_Dict[cpd]['formula']=new_formula
        Update_Formula+=1
        Update_Compounds=1
        
    current_charge = Compounds_Dict[cpd]['charge']
    if(str(new_charge) != current_charge):
        print("Updating charge for "+cpd+" with "+str(new_charge))
        Compounds_Dict[cpd]['charge']=str(new_charge)
        Update_Charge+=1
        Update_Compounds=1

if(Update_Compounds>0):
    print("Updating "+str(Update_Formula)+" compounds with formula and "+str(Update_Charge)+" compounds with charge")
    CompoundsHelper.saveCompounds(Compounds_Dict)

#(inchi_formula,inchi_layers) = InChIs.parse(Structures_Dict[cpd]['InChI'])
#(inchi_formula, notes) = Compounds.mergeFormula(inchi_formula)
#(adjusted_inchi_formula, notes) = InChIs.adjust_protons(inchi_formula, inchi_layers['p'])
