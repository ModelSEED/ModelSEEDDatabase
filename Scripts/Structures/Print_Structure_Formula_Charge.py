#!/usr/bin/env python
import os, sys, re
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

from openbabel import pybel
from rdkit.Chem import AllChem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
resolved_structures=open('Resolved_Structures.txt','w')
unresolved_structures=open('Unresolved_Structures.txt','w')
file_handle_dict=dict()
for source in "KEGG","MetaCyc","Rhea":
    for struct_type in "InChI","SMILE":
        for struct_stage in "Charged","Original":
            file_string="_".join((source,struct_type,struct_stage))
            output_file_name=Structures_Root+source+"/"+struct_type+"_"+struct_stage+"_Formulas_Charges.txt"
            ofh = open(output_file_name,'w')

            input_file_name=Structures_Root+source+"/"+struct_type+"_"+struct_stage+"Strings.txt"
            with open(input_file_name) as ifh:
                for line in ifh.readlines():
                    line = line.strip('\r\n')
                    (external_id,structure)=line.split('\t')

                    mol_rdkit=None
                    mol_openbabel=None
                    try:
                        if(struct_type == 'InChI'):
                            mol_rdkit = AllChem.MolFromInchi(structure)
                            mol_openbabel=pybel.readstring("inchi",structure)
                        elif(struct_type == 'SMILE'):
                            mol_rdkit = AllChem.MolFromSmiles(structure)
                            mol_openbabel=pybel.readstring("smiles",structure)
                    except Exception as e:
                        pass

                    if(mol_rdkit is None):
                        unresolved_structures.write(external_id+"\t"+struct_stage+"\t"+structure+"\tRDKit\n")
                    if(mol_openbabel is None):
                        unresolved_structures.write(external_id+"\t"+struct_stage+"\t"+structure+"\tOpenBabel\n")
                    if(mol_rdkit is None and mol_openbabel is None):
                        continue

                    new_formula=""
                    new_charge=0
                    mol_source=""
                    if(mol_rdkit is not None):
                        mol_source="RDKit"
                        new_formula = AllChem.CalcMolFormula(mol_rdkit)
                        new_charge = AllChem.GetFormalCharge(mol_rdkit)
                        match = re.search('([-+]\d?)$',new_formula)
                        if(match):
                            new_formula = new_formula.replace(match.group(),'')
                    elif(mol_openbabel is not None):
                        mol_source="OpenBabel"
                        new_formula=mol_openbabel.formula
                        new_charge = mol_openbabel.charge
                        match = re.search('([-+]+)$',new_formula)
                        if(match):
                            new_formula = new_formula.replace(match.group(),'')
                    
                    #For SMILES, generic groups are given a '*' character
                    #We're normalizing these as 'R' groups in MSD
                    norm_formula = re.sub('\*','R',new_formula)

                    #normalizing formula my own way, so I can be consistent
                    #these are hill-sorted, and merges molecular fragments
                    norm_formula = Compounds.mergeFormula(norm_formula)[0]

                    resolved_structures.write("\t".join([external_id,struct_stage,structure,norm_formula,str(new_charge),mol_source])+"\n")
                    ofh.write("\t".join((external_id,norm_formula,str(new_charge)))+"\n")