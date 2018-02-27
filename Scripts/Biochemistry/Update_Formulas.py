#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Structures_Dict = CompoundsHelper.loadStructures(["InChI"],["ModelSEED"])

Update_Compounds=0
for cpd in sorted(Compounds_Dict.keys()):
    if(cpd not in Structures_Dict or 'InChI' not in Structures_Dict[cpd]):
        continue

    current_formula = Compounds_Dict[cpd]['formula']

    (inchi_formula,inchi_layers) = InChIs.parse(Structures_Dict[cpd]['InChI'])
    (inchi_formula, notes) = Compounds.mergeFormula(inchi_formula)
    (adjusted_inchi_formula, notes) = InChIs.adjust_protons(inchi_formula, inchi_layers['p'])

    if(adjusted_inchi_formula != current_formula):
        Compounds_Dict[cpd]['formula']=adjusted_inchi_formula
        Update_Compounds+=1

if(Update_Compounds>0):
    print "Saving formula for "+str(Update_Compounds)+" compounds";
    CompoundsHelper.saveCompounds(Compounds_Dict)
