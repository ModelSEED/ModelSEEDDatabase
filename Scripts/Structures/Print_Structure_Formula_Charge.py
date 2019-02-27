#!/usr/bin/env python
import os, sys, re
temp=list();
header=1;

from BiochemPy import Compounds

import pybel
from rdkit.Chem import AllChem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

#Load Structures and Aliases
CompoundsHelper = Compounds()
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChI"],["KEGG","MetaCyc"])

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
file_handle_dict=dict()
for source in "KEGG","MetaCyc":
    for struct_type in "InChI","SMILE":
        for struct_stage in "Charged","Original":
            file_string="_".join((source,struct_type,struct_stage))
            file_name=Structures_Root+source+"/"+struct_type+"_"+struct_stage+"_Formulas_Charges.txt"
            file_handle_dict[file_string]=open(file_name,"w")

unresolved_structures=open('Unresolved_Structures.txt','w')
for struct_type in sorted(Structures_Dict.keys()):
    for external_id in sorted(Structures_Dict[struct_type]):
        source="MetaCyc"
        if(re.search("^[CR]\d{5}$",external_id)):
            source="KEGG"
        for struct_stage in sorted(Structures_Dict[struct_type][external_id].keys()):
            file_string="_".join((source,struct_type,struct_stage))
            for structure in sorted(Structures_Dict[struct_type][external_id][struct_stage].keys()):
                mol=None
                mol_source=""
                try:
                    if(struct_type == 'InChI'):
                        mol = AllChem.MolFromInchi(structure)
                        if(mol is None):
                            mol=pybel.readstring("inchi",structure)
                            if(mol):
                                mol_source="OpenBabel"
                        else:
                            mol_source="RDKit"
                    elif(struct_type == 'SMILE'):
                        mol = AllChem.MolFromSmiles(structure)
                        if(mol==None):
                            mol=pybel.readstring("smiles",structure)
                            if(mol):
                                mol_source="OpenBabel"
                        else:
                            mol_source="RDKit"
                except Exception as e:
                    pass

                if(mol is None):
                    unresolved_structures.write(external_id+"\t"+struct_stage+"\t"+structure+"\n")
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

                #For SMILES, generic groups are given a '*' character
                #We're normalizing these as 'R' groups in MSD
                norm_formula = re.sub('\*','R',new_formula)

                #normalizing formula my own way, so I can be consistent
                #these are hill-sorted, and merges molecular fragments
                norm_formula = Compounds.mergeFormula(norm_formula)[0]

                file_handle_dict[file_string].write("\t".join((external_id,norm_formula,str(new_charge)))+"\n")
