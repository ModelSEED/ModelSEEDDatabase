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
    if(cpd not in Structures_Dict):
        continue

    if('InChI' not in Structures_Dict[cpd]):
        continue

    current_charge = int(Compounds_Dict[cpd]['charge'])

    #Parse out InChI formula and layers
    (inchi_formula,inchi_layers) = InChIs.parse(Structures_Dict[cpd]['InChI'])
    inchi_charge = InChIs.charge(inchi_layers['q'],inchi_layers['p'])

    if(inchi_charge != current_charge):
        Compounds_Dict[cpd]['charge']=str(inchi_charge)
        Update_Compounds+=1

if(Update_Compounds>0):
    print "Saving charge for "+str(Update_Compounds)+" compounds";
    CompoundsHelper.saveCompounds(Compounds_Dict)
