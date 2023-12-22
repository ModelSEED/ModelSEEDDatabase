#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
Structures_Dict = CompoundsHelper.loadStructures(["InChI"],["ModelSEED"])

diff_file = open("Compound_Charge_Differences.txt", 'w')
for cpd in sorted(Compounds_Dict.keys()):
    if(cpd not in Structures_Dict):
        #diff_file.write("Zero structures for "+cpd+"\n")
        continue

    if('InChI' not in Structures_Dict[cpd]):
        #diff_file.write("No InChI structure for "+cpd+"\n")
        continue

    current_charge = float(Compounds_Dict[cpd]['charge'])

    #Parse out InChI formula and layers
    inchi = list(Structures_Dict[cpd]['InChI'].keys())[0]
    (inchi_formula,inchi_layers) = InChIs.parse(inchi)

    inchi_charge = InChIs.charge(inchi_layers['q'],inchi_layers['p'])

    if(inchi_charge != current_charge):
        #Proton-specific (i.e. minor difference)
        if(inchi_layers['q'] == ""):
            diff_file.write("Proton difference for "+cpd+": "+str(current_charge)+" / "+str(inchi_charge)+"\n")
        else:
            diff_file.write("Charge difference for "+cpd+": "+str(current_charge)+" / "+str(inchi_charge)+"\n")
