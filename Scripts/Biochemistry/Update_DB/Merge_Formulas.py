#!/usr/bin/env python

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import os, sys
temp=list();
header=1;

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

Update_Compounds=0
for cpd in sorted(Compounds_Dict.keys()):
    # Skip Light!
    if(cpd == 'cpd11632'):
        continue
    
    old_formula=Compounds_Dict[cpd]["formula"]
    (new_formula, notes) = CompoundsHelper.mergeFormula(old_formula)

    if(notes != ""):
        Compounds_Dict[cpd]["notes"]=notes
        Update_Compounds=1

    if(new_formula != old_formula):
        print("Updating "+cpd+": "+old_formula+" --> "+new_formula)
        Compounds_Dict[cpd]["formula"]=new_formula
        Update_Compounds=1

if(Update_Compounds==1):
    print("Saving compounds")
    CompoundsHelper.saveCompounds(Compounds_Dict)
