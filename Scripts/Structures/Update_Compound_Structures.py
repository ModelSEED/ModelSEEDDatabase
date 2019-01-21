#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

CompoundsHelper = Compounds()
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey"],["ModelSEED"])
Compounds_Dict = CompoundsHelper.loadCompounds()

for cpd in sorted (Compounds_Dict.keys()):
    if(cpd not in Structures_Dict):
        Compounds_Dict[cpd]['inchikey']=""
        Compounds_Dict[cpd]['smiles']=""
    else:
        Compounds_Dict[cpd]['inchikey']=Structures_Dict[cpd].get('InChIKey',"")
        Compounds_Dict[cpd]['smiles']=Structures_Dict[cpd].get('SMILE',"")

print("Saving compounds")
CompoundsHelper.saveCompounds(Compounds_Dict)
