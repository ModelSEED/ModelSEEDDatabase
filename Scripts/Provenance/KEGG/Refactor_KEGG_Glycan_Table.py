#!/usr/bin/env python
import re, sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Compounds

##########################################
# NB: No structures for glycans
##########################################

glycan_file = "data/glycan"
glycan_dict = dict()
glycans_dict = list()
with open(glycan_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_list=re.split('\s+',line)
		data=" ".join(tmp_list[1:])
		if(tmp_list[0] != ""):
			field = tmp_list[0]

		if(field == "ENTRY"):
			glycan_dict = {'id':tmp_list[1],
							'names':[],
							'formula':"null",
							'mass':"10000000",
							'charge':"10000000",
							'inchikey':'',
							'smiles':''}

		elif(field == "NAME"):
			name=data.rstrip(';')
			glycan_dict['names'].append(name)
		elif(field == "COMPOSITION"):
			if(len(glycan_dict['names'])==0):
				name = re.sub('\s','',data)
				glycan_dict['names'].append(name)
		else:
			#print(field)
			pass

		if(field == "///"):

			#print(glycan_dict)
			glycans_dict.append(glycan_dict)

with open("KEGG_glycans.tsv",'w') as kcfh:
	kcfh.write("\t".join(["id","names","formula","mass","charge","inchikey","smiles"])+"\n")
	for gly in glycans_dict:
		kcfh.write("\t".join([gly['id'],"|".join(gly['names']),gly['formula'], \
								gly['mass'],gly['charge'],gly['inchikey'],gly['smiles']])+"\n")