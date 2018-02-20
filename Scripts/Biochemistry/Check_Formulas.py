#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds, InChIs

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

for cpd in sorted(Compounds_Dict.keys()):
    if(Compounds_Dict[cpd]['inchikey'] == 'null'):
        continue

    inchi_formula = InChIs.parse(Compounds_Dict[cpd]['inchikey'],True)

    if(inchi_formula[0][0] != Compounds_Dict[cpd]['formula']):
        print cpd, inchi_formula, Compounds_Dict[cpd]['formula']
        break



#Update_Compounds=0
#for cpd in sorted(Compounds_Dict.keys()):
#    old_formula=Compounds_Dict[cpd]["formula"]
#    (new_formula, notes) = CompoundsHelper.mergeFormula(old_formula)

#    if(notes != ""):
#        Compounds_Dict[cpd]["notes"]=notes
#        Update_Compounds=1

#    if(new_formula != old_formula):
#        print "Updating "+cpd+": "+old_formula+" --> "+new_formula
#        Compounds_Dict[cpd]["formula"]=new_formula
#        Update_Compounds=1

#if(Update_Compounds==1):
#    print "Saving componds";
#    CompoundsHelper.saveCompounds(Compounds_Dict)
