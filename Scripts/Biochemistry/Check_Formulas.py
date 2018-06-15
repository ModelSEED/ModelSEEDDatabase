#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Structures_Dict = CompoundsHelper.loadStructures(["InChI"],["ModelSEED"])

for cpd in sorted(Compounds_Dict.keys()):
    if(Compounds_Dict[cpd]['inchikey'] == '' or cpd not in Structures_Dict):
        continue

    (inchi_formula,inchi_layers) = InChIs.parse(Structures_Dict[cpd]['InChI'])
    merged_inchi_formula = CompoundsHelper.mergeFormula(inchi_formula)[0]
    adjusted_inchi_formula = (InChIs.adjust_protons(merged_inchi_formula,inchi_layers['p']))[0]

    if(adjusted_inchi_formula != Compounds_Dict[cpd]['formula']):

         adjusted_inchi_atoms_dict = CompoundsHelper.parseFormula(adjusted_inchi_formula)
         if('H' in adjusted_inchi_atoms_dict):
             del(adjusted_inchi_atoms_dict['H'])
         adjusted_inchi_protonfree_formula = CompoundsHelper.buildFormula(adjusted_inchi_atoms_dict)

         original_formula_atoms_dict = CompoundsHelper.parseFormula(Compounds_Dict[cpd]['formula'])
         if('H' in original_formula_atoms_dict):
             del(original_formula_atoms_dict['H'])
         original_formula_protonfree_formula = CompoundsHelper.buildFormula(original_formula_atoms_dict)

         if(adjusted_inchi_protonfree_formula != original_formula_protonfree_formula):
             print cpd,original_formula_protonfree_formula,adjusted_inchi_protonfree_formula
#             break
