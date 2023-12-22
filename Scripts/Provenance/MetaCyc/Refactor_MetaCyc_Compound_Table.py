#!/usr/bin/env python
import re, sys
sys.path.append('../../../Libs/Python/')
from BiochemPy import Compounds

##########################################
# Load inchikey structures
MetaCyc_structures_dict = dict()

MetaCyc_File = "../../../Biochemistry/Structures/MetaCyc/InChIKey_ChargedStrings.txt"
with open(MetaCyc_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in MetaCyc_structures_dict):
			MetaCyc_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		MetaCyc_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

MetaCyc_File = "../../../Biochemistry/Structures/MetaCyc/InChIKey_OriginalStrings.txt"
with open(MetaCyc_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in MetaCyc_structures_dict):
			MetaCyc_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		
		if(MetaCyc_structures_dict[tmp_list[0]]['inchikey'] == ''):
			MetaCyc_structures_dict[tmp_list[0]]['inchikey']=tmp_list[1]

##########################################
# Load smile structures
MetaCyc_File = "../../../Biochemistry/Structures/MetaCyc/SMILE_ChargedStrings.txt"
with open(MetaCyc_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in MetaCyc_structures_dict):
			MetaCyc_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		MetaCyc_structures_dict[tmp_list[0]]['smiles']=tmp_list[1]

MetaCyc_File = "../../../Biochemistry/Structures/MetaCyc/SMILE_OriginalStrings.txt"
with open(MetaCyc_File) as fh:
	for line in fh.readlines():
		line = line.strip('\r\n')
		tmp_list = line.split('\t')
		if(tmp_list[0] not in MetaCyc_structures_dict):
			MetaCyc_structures_dict[tmp_list[0]]={'inchikey':'','smiles':''}
		elif(MetaCyc_structures_dict[tmp_list[0]] == ''):
			MetaCyc_structures_dict[tmp_list[0]]['smiles']=tmp_list[1]

def convert_html_name(name):

	if('&' in name):
		
		#greek letters
		for greek in ['alpha','beta','gamma','delta','epsilon','omega', \
						'mu','nu','kappa','chi','zeta','psi','pi','phi', \
							'tau','iota', 'theta','sigma','lambda','xi']:
			name = re.sub('&'+greek+';',greek,name)
			Greek = greek[0].upper()+greek[1:]
			name = re.sub('&'+Greek+';',greek,name)

		for html in ['i','b','a','em','small','sup','sub','span','href']:
			name = re.sub('</?'+html+'/?>','',name)
			HTML = html.upper()
			name = re.sub('</?'+HTML+'/?>','',name)

		#arrows
		name = re.sub('&rarr;','->',name)
		name = re.sub('&RARR;','->',name)
		name = re.sub('&harr;','<->',name)

		#dashes
		name = re.sub('&ndash;','-',name)
		name = re.sub('&mdash;','-',name)

		#middle dot
		name = re.sub('&middot;','-',name)

		#plus/minus
		name = re.sub('&plusmn;','+/-',name)
		
		#ampersand?!
		name = re.sub('&amp;','',name)

		if('&' in name or ';' in name):
			print("Warning: HTML characters in "+name)

		name = re.sub(';',',',name)

	return name

"""
    $bit =~ s/\&plusmn;/\(+\/-\)/g;
    $bit =~ s/\&prime;?/\'/g;
    $bit =~ s/\&deg;//g;
    $bit =~ s/\&acirc;?/y/g;
    $bit =~ s/\&mdash;/\-/g;
    $bit =~ s/\&\#8596;/\<\-\>/g;
    $bit =~ s/\&\#8217;?/\'/g;
    $bit =~ s/\&\#8242;/\'/g;
    $bit =~ s/\&\#8222;/\'/g;
"""

compound_dict = dict()
compounds_list = list()
with open('data/compounds.dat') as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')

		# Skip commented lines
		if(line.startswith('#')):
			continue

		tmp_list = line.split(' - ')
		field = tmp_list[0]
		data="-".join(tmp_list[1:])
		data=data.strip(' ')

		if(field == "UNIQUE-ID"):

			compound_dict = {'id':data,
							'names':[],
							'formula':"null",
							'mass':"10000000",
							'charge':"10000000",
							'inchikey':'',
							'smiles':''}

		elif(field == "CHEMICAL-FORMULA"):

			# Strip parentheses
			formula=data.strip('(')
			formula=formula.strip(')')

			# Split element and stoichiometry
			[element,stoichiometry]=formula.split(' ')

			# Reduce second letter of element to lowercase
			if(len(element)==2):
				element = element[0]+element[1].lower()

			# Warn about weird element
			if(len(element)>2):
				print("Warning, un-natural element "+element+" for compound "+compound_dict['id'])
			else:

				formula = element
				if(stoichiometry != '1'):
					formula+=stoichiometry
				if(compound_dict['formula'] == 'null'):
					compound_dict['formula']=formula
				else:
					compound_dict['formula']+=formula

		elif(field == "SYNONYMS" or field.endswith('NAME')):
			if(field == "N-NAME" or field == "N+1-NAME" or field == "N-1-NAME"):
				continue

			name = convert_html_name(data)
			if(name not in compound_dict['names']):
				compound_dict['names'].append(name)

		elif(field == "CHARGE"):
			compound_dict['charge']=data

		elif(field == "MOLECULAR-WEIGHT"):
			mass = re.sub('\.$','',data)
			compound_dict['mass']=mass

		elif(field == "SMILES"):
			compound_dict['smiles']=data
		elif(field == 'INCHI'):
			pass
		elif(field.startswith("//")):
			compound_dict['formula'] = Compounds.mergeFormula(compound_dict['formula'])[0]

			#Adjust for generic compounds (and not Radon)
			if('R' in compound_dict['formula'] and 'Rn' not in compound_dict['formula']):
				compound_dict['mass']="0"
				compound_dict['charge']="0"

			#Modifictions for photon
			if(compound_dict['id'] == "Light" or compound_dict['id'] == 'UV-Light'):
				compound_dict['mass']="0"
				compound_dict['charge']="0"
	
			#Populate with structures
			if(compound_dict['id'] in MetaCyc_structures_dict):
				compound_dict['inchikey']=MetaCyc_structures_dict[compound_dict['id']]['inchikey']
				compound_dict['smiles']=MetaCyc_structures_dict[compound_dict['id']]['smiles']

			compounds_list.append(compound_dict)

with open("MetaCyc_compounds.tsv",'w') as mcfh:
	mcfh.write("\t".join(["id","names","formula","mass","charge","inchikey","smiles"])+"\n")
	for cpd in compounds_list:
		mcfh.write("\t".join([cpd['id'],"|".join(cpd['names']),cpd['formula'], \
								cpd['mass'],cpd['charge'],cpd['inchikey'],cpd['smiles']])+"\n")
