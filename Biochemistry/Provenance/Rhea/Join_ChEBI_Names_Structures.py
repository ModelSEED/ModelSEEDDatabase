#!/usr/bin/env python
##########################################
# Load names
ChEBI_names_dict = dict()
header=1
with open('ChEBI_rdf_names.tsv') as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(header==1):
			header=0
			continue
		chebi_id = tmp_list[0]
		if(chebi_id not in ChEBI_names_dict):
			ChEBI_names_dict[chebi_id]=list()
		ChEBI_names_dict[chebi_id].append(tmp_list[1])

with open('names.tsv') as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[1] not in ChEBI_names_dict):
			continue
		if(tmp_list[4] not in ChEBI_names_dict[tmp_list[1]]):
			ChEBI_names_dict[tmp_list[1]].append(tmp_list[4])

##########################################
# Load inchikey structures
ChEBI_structures_dict = dict()

ChEBI_File = "../../../Structures/ChEBI/InChIKey_ChargedStrings.txt"
with open(ChEBI_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in ChEBI_structures_dict):
			ChEBI_structures_dict[tmp_list[0]]={'inchikey':'','smile':''}
		ChEBI_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

ChEBI_File = "../../../Structures/ChEBI/InChIKey_OriginalStrings.txt"
with open(ChEBI_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in ChEBI_structures_dict):
			ChEBI_structures_dict[tmp_list[0]]={'inchikey':'','smile':''}
		
		if(ChEBI_structures_dict[tmp_list[0]]['inchikey'] == ''):
			ChEBI_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

##########################################
# Load smile structures
ChEBI_File = "../../../Structures/ChEBI/SMILE_ChargedStrings.txt"
with open(ChEBI_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in ChEBI_structures_dict):
			ChEBI_structures_dict[tmp_list[0]]={'inchikey':'','smile':''}
		ChEBI_structures_dict[tmp_list[0]]['smile']=tmp_list[1]

ChEBI_File = "../../../Structures/ChEBI/SMILE_OriginalStrings.txt"
with open(ChEBI_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in ChEBI_structures_dict):
			ChEBI_structures_dict[tmp_list[0]]={'inchikey':'','smile':''}
		elif(ChEBI_structures_dict[tmp_list[0]] == ''):
			ChEBI_structures_dict[tmp_list[0]]['smile']=tmp_list[1]

TSV_File = "ChEBI_ID_Name_Structure.tsv"
with open(TSV_File, 'w') as ofh:
	ofh.write('\t'.join(['ID','Names','InChIKey','SMILES'])+'\n')
	for chebi_id in sorted(ChEBI_names_dict):
		line_list = [chebi_id]
		line_list.append('|'.join(ChEBI_names_dict[chebi_id]))
		if(chebi_id in ChEBI_structures_dict):
			line_list.append(ChEBI_structures_dict[chebi_id]['inchikey'])
			line_list.append(ChEBI_structures_dict[chebi_id]['smile'])
		else:
			line_list.append('')
			line_list.append('')

		line='\t'.join(line_list)
		ofh.write(line+'\n')