#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Structures_Dict = CompoundsHelper.loadStructures(["InChI"],["ModelSEED"])

diff_file = open("Compound_Formula_Differences.txt", 'w')
for cpd in sorted(Compounds_Dict.keys()):
    if(cpd not in Structures_Dict):
        diff_file.write("Zero structures for "+cpd+"\n")
        continue

    if('InChI' not in Structures_Dict[cpd]):
        diff_file.write("No InChI structure for "+cpd+"\n")
        continue

    current_formula = Compounds_Dict[cpd]['formula']

    #Parse out InChI formula
    (inchi_formula,inchi_layers) = InChIs.parse(Structures_Dict[cpd]['InChI'])

    #Make sure formula is merged appropriately before applying proton adjustment
    (inchi_formula, notes) = Compounds.mergeFormula(inchi_formula)
    if(notes != ""):
        diff_file.write("Notes from merging InChI formula for "+cpd+": "+notes+"\n")

    #Make adjustments based to protonation state of InChI structure
    (adjusted_inchi_formula, notes) = InChIs.adjust_protons(inchi_formula, inchi_layers['p'])
    if(notes != ""):
        diff_file.write("Notes from adjusting protons for "+cpd+": "+notes+"\n")

    if(adjusted_inchi_formula != current_formula):
        if(current_formula == "null"):
            diff_file.write("New formula for "+cpd+": "+adjusted_inchi_formula+"\n")

        else:

            #Check to see if difference in update is only in protons
            old_atoms = Compounds.parseFormula(current_formula)
            new_atoms = Compounds.parseFormula(adjusted_inchi_formula)

            missing_atoms=set()
            different_atoms=dict()
            for atom in old_atoms.keys():
                if(atom not in new_atoms):
                    missing_atoms.update(atom)
                elif(old_atoms[atom] != new_atoms[atom]):
                    different_atoms[atom]=old_atoms[atom]-new_atoms[atom]

            for atom in new_atoms.keys():
                if(atom not in old_atoms):
                    missing_atoms.add(atom)

            #Proton-specific (i.e. minor difference)
            if( (len(missing_atoms)==1 and 'H' in missing_atoms) or (len(different_atoms)==1 and "H" in different_atoms)):
                diff_file.write("Proton difference for "+cpd+": "+str(missing_atoms)+"/"+str(different_atoms)+"\n")
                continue

            if(len(missing_atoms)>0):
                diff_file.write("Missing atoms for "+cpd+": "+str(missing_atoms)+"\n")
            if(len(different_atoms.keys())>0):
                diff_file.write("Differing atoms for "+cpd+": "+str(different_atoms)+"\n")
