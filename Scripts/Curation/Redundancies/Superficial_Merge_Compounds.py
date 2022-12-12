#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../../Libs/Python')
from BiochemPy import Reactions, Compounds

arguments = list(sys.argv)
#Remove script name
arguments = arguments[1:]
if(len(arguments) != 2 or 'cpd' not in arguments[0] or 'cpd' not in arguments[1]):
    print("Error: script must be initiated with the identifiers of the two compounds to merge")
    sys.exit()

primary_cpd=arguments[0]
merging_cpd=arguments[1]
arguments=sorted(arguments)

if(primary_cpd != arguments[0]):
    print("Error: compound identifiers must be used in order")
    print("\tThe first compound identifier (in order) should be the one that is retained: "+arguments[0])
    sys.exit()

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

if(compounds_dict[primary_cpd]['is_obsolete']==1):
    compounds_dict[primary_cpd]["is_obsolete"]=0

if(merging_cpd not in compounds_dict[primary_cpd]['linked_compound']):
    lnkd_cpds = compounds_dict[primary_cpd]['linked_compound'].split(';')
    if('null' in lnkd_cpds):
        lnkd_cpds.remove('null')
    lnkd_cpds.append(merging_cpd)
    compounds_dict[primary_cpd]['linked_compound']=";".join(sorted(lnkd_cpds))

if(compounds_dict[merging_cpd]['is_obsolete']==0):
    compounds_dict[merging_cpd]['is_obsolete']=1

if(primary_cpd not in compounds_dict[merging_cpd]['linked_compound']):
    lnkd_cpds = compounds_dict[merging_cpd]['linked_compound'].split(';')
    if('null' in lnkd_cpds):
        lnkd_cpds.remove('null')
    lnkd_cpds.append(primary_cpd)
    compounds_dict[merging_cpd]['linked_compound']=";".join(sorted(lnkd_cpds))

compounds_helper.saveCompounds(compounds_dict)

print("You must run these commands afterwards:")
print("\t../../Biochemistry/Maintain/Update_Obsolete_Compounds_in_Reactions.py")
print("\t../../Biochemistry/Refresh/Refresh_Aliases.sh")
print("\t../../Biochemistry/Refresh/Refresh_Reactions.sh")
print("\t../../Biochemistry/Update/Remove_Newly_Obsolescent_Reactions.py")
print("\t../../Biochemistry/Update/Remove_Newly_Obsolescent_Compounds.py")
print("\t../../Biochemistry/Reprint_Biochemistry.py")
print("\t../Structures/List_ModelSEED_Structures.py")
print("\t../Structures/Update_Structure_Formula_Charge.py")
