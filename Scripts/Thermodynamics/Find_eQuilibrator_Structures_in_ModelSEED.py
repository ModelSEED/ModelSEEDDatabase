#!/usr/bin/env python

if __name__ == "__main__":
    # Argument guard -- see "The argument guard" in Scripts/README.md.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


import os

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."

ms_structures_list=list()
with open(ROOT + '/Biochemistry/Structures/All_ModelSEED_Structures.txt') as ms_fh:
    for line in ms_fh.readlines():
        line=line.strip('\r\n')
        tmp_list=line.split('\t')

        if(tmp_list[1] != "InChIKey"):
           continue

        if(tmp_list[7] not in ms_structures_list):
            ms_structures_list.append(tmp_list[7])

#Contrived from equilibrator's cache, see Biochemistry/Structures/MetaNetX/README.md
#08/31/23
#Reads eq_cpds.tsv from the current working directory — drop it next
#to wherever you run this script (the MetaNetX README walks through
#the SQLite export step).
ofh = open(ROOT + '/Biochemistry/Structures/MetaNetX/Structures_in_ModelSEED_and_eQuilibrator.txt','w')
with open('eq_cpds.tsv') as eq_fh:
    for line in eq_fh.readlines():
        line=line.strip('\r\n')
        tmp_list = line.split('\t')

        if(tmp_list[0] == '' or tmp_list[1] == ''):
            continue

        if(tmp_list[1] in ms_structures_list):
            ofh.write('\t'.join(tmp_list[0:2])+'\n')
