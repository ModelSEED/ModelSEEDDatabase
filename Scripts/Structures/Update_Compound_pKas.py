#!/usr/bin/env python
import os,sys
from BiochemPy import Compounds

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"

# Load pKas and pKbs
cpd_pKab_dict=dict()
for DB in ["KEGG","MetaCyc"]:
    with open(Structures_Root+DB+'/pKa_Strings.txt') as fh:
        for line in fh.readlines():
            line=line.strip()
            array=line.split('\t')
            if(array[0] not in cpd_pKab_dict):
                cpd_pKab_dict[array[0]]={array[1]:array[2]}
            else:
                cpd_pKab_dict[array[0]][array[1]]=array[2]

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
structures_dict = compounds_helper.loadStructures(["SMILE","InChI","InChIKey"],["ModelSEED"])
aliases_dict = compounds_helper.loadMSAliases()

# We're removing all pKa and pKb before loading new ones
for cpd in compounds_dict:
    compounds_dict[cpd]['pka']=""
    compounds_dict[cpd]['pkb']=""

# We're only loading pKa/pKb for compounds that have an accepted unique structure in ModelSEED
for cpd in structures_dict:
    found=False
    for DB in ["KEGG","MetaCyc"]:
        if(found is True or DB not in aliases_dict[cpd]):
            continue

        for alias in aliases_dict[cpd][DB]:
            if(alias in cpd_pKab_dict):
                print(cpd,alias,cpd_pKab_dict[alias])
                if('pKa' in cpd_pKab_dict[alias]):
                    print(cpd,alias)
                    compounds_dict[cpd]['pka']=cpd_pKab_dict[alias]['pKa']
                else:
                    compounds_dict[cpd]['pka']=""

                if('pKb' in cpd_pKab_dict[alias]):
                    compounds_dict[cpd]['pkb']=cpd_pKab_dict[alias]['pKb']
                else:
                    compounds_dict[cpd]['pkb']=""

                # All structures for the same compound should be the same,
                # and so only need to process once
                found=True
                break

print("Saving compounds")
compounds_helper.saveCompounds(compounds_dict)
