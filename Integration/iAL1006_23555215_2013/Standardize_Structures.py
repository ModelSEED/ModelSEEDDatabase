#!/usr/bin/env python
import os, sys, re
temp=list();
header=1;

import pybel
from rdkit.Chem import AllChem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

if(len(sys.argv)==1):
    print("Must pass compounds file with InChI column")
    sys.exit()

Biochem_File = sys.argv[1]
import os.path
if(os.path.exists(Biochem_File) is not True):
    print("No File!")
    sys.exit()

File_Lines=list()
with open(Biochem_File) as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        array=line.split('\t')
        if(array[2] != "InChI" and array[2] != ""):

            (mol_source,inchi_string)=("","")
            
            try:
                mol = AllChem.MolFromInchi(array[2])
                if(mol is None):
                    mol=pybel.readstring("inchi",array[2])
                    if(mol):
                        mol_source="OpenBabel"
                else:
                    mol_source="RDKit"
            except Exception as e:
                print(e)
                pass

            if(mol_source=="RDKit"):
                (inchi_string,AuxInfo) = AllChem.MolToInchiAndAuxInfo(mol)
            elif(mol_source=="OpenBabel"):
                inchi_string=mol.write(format="inchi").strip()

            if(inchi_string==""):
                print("Failed to convert InChI string for "+array[0])
            else:
                array[2]=inchi_string
        
        File_Lines.append('\t'.join(array))
fh.close()

with(open(Biochem_File,'w')) as fh:
    for line in File_Lines:
        fh.write(line+'\n')
fh.close()
