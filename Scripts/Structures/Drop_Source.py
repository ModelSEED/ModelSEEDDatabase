#!/usr/bin/env python
Sources = ['ChEBI','Rhea']
IDs_to_remove = list()
All_lines = list()
with open('../../Biochemistry/Structures/All_ModelSEED_Structures.txt') as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[4] not in Sources):
			All_lines.append(line)
		else:
			IDs_to_remove.append(tmp_list[3])

with open('../../Biochemistry/Structures/All_ModelSEED_Structures.txt','w') as fh:
	fh.write('\n'.join(All_lines)+'\n')

All_lines = list()
with open('../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt') as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		ids_list = tmp_list[2].split(';')
		new_ids_list = list()
		for id in ids_list:
			if(id == ''):
				continue

			if(id not in IDs_to_remove):
				new_ids_list.append(id)
		tmp_list[2] = ";".join(new_ids_list)
		line = "\t".join(tmp_list)
		if(len(new_ids_list)>0):
			All_lines.append(line)

with open('../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt', 'w') as fh:
	fh.write('\n'.join(All_lines)+'\n')