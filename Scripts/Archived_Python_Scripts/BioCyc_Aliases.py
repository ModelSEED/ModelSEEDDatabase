#!/usr/bin/env python
import os
import sys
import json
from BiochemPy import Compounds, Reactions

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["MetaCyc","PlantCyc"])

for cpd in MS_Aliases_Dict:
    if('PlantCyc' not in MS_Aliases_Dict[cpd]):
        MS_Aliases_Dict[cpd]['PlantCyc']=[]
    if('MetaCyc' not in MS_Aliases_Dict[cpd]):
        MS_Aliases_Dict[cpd]['MetaCyc']=[]

    print("\t".join([cpd,"|".join(MS_Aliases_Dict[cpd]['PlantCyc']),"|".join(MS_Aliases_Dict[cpd]["MetaCyc"])]))

#Load Reactions
ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
MS_Aliases_Dict =  ReactionsHelper.loadMSAliases(["MetaCyc","PlantCyc"])

for rxn in MS_Aliases_Dict:
    if('PlantCyc' not in MS_Aliases_Dict[rxn]):
        MS_Aliases_Dict[rxn]['PlantCyc']=[]
    if('MetaCyc' not in MS_Aliases_Dict[rxn]):
        MS_Aliases_Dict[rxn]['MetaCyc']=[]

    print("\t".join([rxn,"|".join(MS_Aliases_Dict[rxn]['PlantCyc']),"|".join(MS_Aliases_Dict[rxn]["MetaCyc"])]))

