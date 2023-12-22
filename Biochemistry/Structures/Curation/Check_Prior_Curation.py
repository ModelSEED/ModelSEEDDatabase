#!/usr/bin/env python
import glob

#Load Curated Structures
ignored_structures=dict()
for file in glob.glob("*.txt"):
    with open(file) as ignore_file:
        for line in ignore_file.readlines():
            array=line.split('\t')
            ignored_structures[array[0]]=1

file = "../All_ModelSEED_Structures.txt"
ignored_cpd_inchikey=dict()
with open(file) as afh:
	for line in afh.readlines():
		line=line.strip('\r\n')
		tmp_list=line.split('\t')
		if(tmp_list[3] in ignored_structures):
			if(tmp_list[1] == "InChIKey" and tmp_list[2] == "Charged"):
				ignored_cpd_inchikey[tmp_list[0]]={'ick':tmp_list[7],'id':tmp_list[3]}
	
file = "../Unique_ModelSEED_Structures.txt"
with open(file) as afh:
	for line in afh.readlines():
		line=line.strip('\r\n')
		tmp_list=line.split('\t')
		if(tmp_list[0] in ignored_cpd_inchikey and \
			tmp_list[5] in ignored_cpd_inchikey[tmp_list[0]]['ick']):
			print(ignored_cpd_inchikey[tmp_list[0]]['id'])