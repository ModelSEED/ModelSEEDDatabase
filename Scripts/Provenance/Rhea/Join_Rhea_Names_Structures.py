#!/usr/bin/env python

if __name__ == "__main__":
    # Validate arguments BEFORE importing anything or touching the database.
    # These scripts mutate the database, and without this an unknown flag or a
    # mistyped mode was silently ignored and the script ran with its defaults:
    # asking Estimate_Reaction_Reversibility.py for --help rewrote 122 files.
    # Placed above the imports so --help works even where a dependency is
    # missing from the path.
    import argparse as _argparse
    _argparse.ArgumentParser(
        description=__doc__,
        formatter_class=_argparse.RawDescriptionHelpFormatter).parse_args()


##########################################
# Load names
Rhea_names_dict = dict()
header=1
with open('Rhea_rdf_names.tsv') as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(header==1):
			header=0
			continue
		if(tmp_list[0] not in Rhea_names_dict):
			Rhea_names_dict[tmp_list[0]]=list()
		Rhea_names_dict[tmp_list[0]].append(tmp_list[1])

##########################################
# Load inchikey structures
Rhea_structures_dict = dict()

Rhea_File = "../../../Structures/Rhea/InChIKey_ChargedStrings.txt"
with open(Rhea_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		rhea_id = ':'.join(tmp_list[0].split('_'))
		if(rhea_id not in Rhea_structures_dict):
			Rhea_structures_dict[rhea_id]={'inchikey':'','smile':''}
		Rhea_structures_dict[rhea_id]['inchikey']=tmp_list[1]

Rhea_File = "../../../Structures/Rhea/InChIKey_OriginalStrings.txt"
with open(Rhea_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		rhea_id = ':'.join(tmp_list[0].split('_'))
		if(rhea_id not in Rhea_structures_dict):
			Rhea_structures_dict[rhea_id]={'inchikey':'','smile':''}
		
		if(Rhea_structures_dict[rhea_id]['inchikey'] == ''):
			Rhea_structures_dict[rhea_id]['inchikey']=tmp_list[1]

##########################################
# Load smile structures
Rhea_File = "../../../Structures/Rhea/SMILE_ChargedStrings.txt"
with open(Rhea_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		rhea_id = ':'.join(tmp_list[0].split('_'))
		if(rhea_id not in Rhea_structures_dict):
			Rhea_structures_dict[rhea_id]={'inchikey':'','smile':''}
		Rhea_structures_dict[rhea_id]['smile']=tmp_list[1]

Rhea_File = "../../../Structures/Rhea/SMILE_OriginalStrings.txt"
with open(Rhea_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		rhea_id = ':'.join(tmp_list[0].split('_'))
		if(rhea_id not in Rhea_structures_dict):
			Rhea_structures_dict[rhea_id]={'inchikey':'','smile':''}
		elif(Rhea_structures_dict[rhea_id] == ''):
			Rhea_structures_dict[rhea_id]['smile']=tmp_list[1]

TSV_File = "Rhea_ID_Name_Structure.tsv"
with open(TSV_File, 'w') as ofh:
	ofh.write('\t'.join(['ID','Names','InChIKey','SMILES'])+'\n')
	for rhea_id in sorted(Rhea_names_dict):
		line_list = [rhea_id]
		line_list.append('|'.join(Rhea_names_dict[rhea_id]))
		if(rhea_id in Rhea_structures_dict):
			line_list.append(Rhea_structures_dict[rhea_id]['inchikey'])
			line_list.append(Rhea_structures_dict[rhea_id]['smile'])
		else:
			line_list.append('')
			line_list.append('')

		line='\t'.join(line_list)
		ofh.write(line+'\n')