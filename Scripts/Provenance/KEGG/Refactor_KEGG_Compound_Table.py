#!/usr/bin/env python
import re, sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Compounds

##########################################
# Load inchikey structures
KEGG_structures_dict = dict()

KEGG_File = "../../../Biochemistry/Structures/KEGG/InChIKey_ChargedStrings.txt"
with open(KEGG_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in KEGG_structures_dict):
			KEGG_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		KEGG_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

KEGG_File = "../../../Biochemistry/Structures/KEGG/InChIKey_OriginalStrings.txt"
with open(KEGG_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in KEGG_structures_dict):
			KEGG_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		
		if(KEGG_structures_dict[tmp_list[0]]['inchikey'] == ''):
			KEGG_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

##########################################
# Load smile structures
KEGG_File = "../../../Biochemistry/Structures/KEGG/SMILE_ChargedStrings.txt"
with open(KEGG_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in KEGG_structures_dict):
			KEGG_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		KEGG_structures_dict[tmp_list[0]]['smiles']=tmp_list[1]

KEGG_File = "../../../Biochemistry/Structures/KEGG/SMILE_OriginalStrings.txt"
with open(KEGG_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in KEGG_structures_dict):
			KEGG_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		elif(KEGG_structures_dict[tmp_list[0]] == ''):
			KEGG_structures_dict[tmp_list[0]]['smiles']=tmp_list[1]

compound_file = "data/compound"
compound_dict = dict()
compounds_dict = list()
with open(compound_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_list=re.split('\s+',line)
		data=" ".join(tmp_list[1:])
		if(tmp_list[0] != ""):
			field = tmp_list[0]

		if(field == "ENTRY"):
			compound_dict = {'id':tmp_list[1],
							'names':[],
							'formula':"null",
							'mass':"10000000",
							'charge':"10000000",
							'inchikey':'',
							'smiles':''}

		elif(field == "NAME"):
			data=data.rstrip(';')
			compound_dict['names'].append(data)
		elif(field == "FORMULA"):
			formula = Compounds.mergeFormula(data)[0]
			compound_dict['formula']=formula
		elif(field == "EXACT_MASS"):
			compound_dict['mass']=data
		else:
			#print(field)
			pass

		if(field == "///"):

			#Adjust for generic compounds (and not Radon)
			if('R' in formula and 'Rn' not in formula):
				compound_dict['charge']="0"
				compound_dict['charge']="0"

			#Modifictions for electron and photon
			if(compound_dict['id'] == "C05359"):
				compound_dict['mass']="0"
				compound_dict['charge']="-1"
	
			if(compound_dict['id'] == "C00205"):
				compound_dict['mass']="0"
				compound_dict['charge']="0"

			if(compound_dict['id'] in KEGG_structures_dict):
				compound_dict['inchikey']=KEGG_structures_dict[compound_dict['id']]['inchikey']
				compound_dict['smiles']=KEGG_structures_dict[compound_dict['id']]['smiles']

			#print(compound_dict)
			compounds_dict.append(compound_dict)

with open("KEGG_compounds.tsv",'w') as kcfh:
	kcfh.write("\t".join(["id","names","formula","mass","charge","inchikey","smiles"])+"\n")
	for cpd in compounds_dict:
		kcfh.write("\t".join([cpd['id'],"|".join(cpd['names']),cpd['formula'], \
								cpd['mass'],cpd['charge'],cpd['inchikey'],cpd['smiles']])+"\n")
