#!/usr/bin/env python
import re

def convert_html_name(name):

	if('&' in name):
		
		#greek letters
		for greek in ['alpha','beta','gamma','delta','epsilon','omega', \
						'mu','nu','kappa','chi','zeta','psi','pi','phi', \
							'tau','iota', 'theta','sigma','lambda','xi']:
			name = re.sub('&'+greek+';',greek,name)
			Greek = greek[0].upper()+greek[1:]
			name = re.sub('&'+Greek+';',greek,name)

		#arrows
		name = re.sub('&rarr;','->',name)
		name = re.sub('&RARR;','->',name)
		name = re.sub('&harr;','<->',name)

		#dashes
		name = re.sub('&ndash;','-',name)
		name = re.sub('&mdash;','-',name)

		#middle dot
		name = re.sub('&middot;','-',name)

		#prime
		name = re.sub('&prime;','\'',name)

		#plus/minus
		name = re.sub('&plusmn;','+/-',name)
		
		#ampersand?!
		name = re.sub('&amp;','',name)

	for html in ['i','b','a','em','small','sup','sub','span','href']:
		name = re.sub('</?'+html+'/?>','',name)
		HTML = html.upper()
		name = re.sub('</?'+HTML+'/?>','',name)

	name = re.sub(';',',',name)

	if('&' in name or ';' in name):
		print("Warning: HTML characters in "+name)

	return name

types_dict=dict()
field_id=None
with open('data/classes.dat',errors='ignore') as rfh:
	for line in rfh.readlines():
		line=line.strip('\r\n')

		# Skip commented lines
		if(line.startswith('#')):
			continue
		
		tmp_list = line.split(' - ')
		field = tmp_list[0]
		data="-".join(tmp_list[1:])
		data=data.strip(' ')
		
		if(field == "UNIQUE-ID"):
			field_id = data
			if(field_id not in types_dict):
				types_dict[field_id]=[]

		elif(field == "TYPES"):

			typed = data
			if(typed not in types_dict[field_id]):
				types_dict[field_id].append(typed)


reaction_list=list()
with open('data/reactions.dat',errors='ignore') as rfh:
	for line in rfh.readlines():
		line=line.strip('\r\n')

		# Skip commented lines
		if(line.startswith('#')):
			continue
		
		tmp_list = line.split(' - ')
		field = tmp_list[0]
		data="-".join(tmp_list[1:])
		data=data.strip(' ')
		
		if(field == "UNIQUE-ID"):
			if(data not in reaction_list):
				reaction_list.append(data)

pathways_dict = dict()
pathway = None
with open('data/pathways.dat',errors='ignore') as pfh:
	for line in pfh.readlines():
		line=line.strip('\r\n')

		# Skip commented lines
		if(line.startswith('#')):
			continue

		tmp_list = line.split(' - ')
		field = tmp_list[0]
		data="-".join(tmp_list[1:])
		data=data.strip(' ')

		if(field == "UNIQUE-ID"):
			pathway = data
			pathways_dict[pathway] = {'name':[],
						  'reactions':[],
						  'pathways':[],
						  'parent':[]}

		elif(field == "COMMON-NAME"):

			name = convert_html_name(data)
			# if(name != data):
			#	print(pathway_dict['id'],data,name)
			pathways_dict[pathway]['name']=name

		elif(field == "REACTION-LIST"):
			if(data not in reaction_list):
				pathways_dict[pathway]['pathways'].append(data)
			else:
				pathways_dict[pathway]['reactions'].append(data)

		elif(field == "TYPES"):

			if(data not in pathways_dict[pathway]['parent']):
				pathways_dict[pathway]['parent'].append(data)

		elif(field.startswith('//')):
			pass

def add_type_as_pwy(child):
	if(child not in pathways_dict):
		pathways_dict[child]={'name':"",
				      'reactions':[],
				      'pathways':[],
				      'parent':[]}
		
	for parent in types_dict[child]:
		if(parent == "Pathways"):
			return

		if(parent not in pathways_dict[child]['parent']):
			pathways_dict[child]['parent'].append(parent)

		if(parent not in pathways_dict):
			add_type_as_pwy(parent)

pathway_types = dict()
pathways_list = list(pathways_dict.keys())
for pwy in pathways_list:
	for typed in pathways_dict[pwy]['parent']:
		add_type_as_pwy(typed)

with open("MetaCyc_pathways.tsv",'w') as mpfh:
	mpfh.write("\t".join(["id","name","reactions","parent"])+"\n")
	for pwy in pathways_dict:
		pwy_data = pathways_dict[pwy]
		for subpwy in pwy_data['pathways']:
			
			if(subpwy not in pathways_dict):
				print("Warning! pathway not found: "+subpwy)
				continue

			for rxn in pathways_dict[subpwy]['reactions']:
				# only occurs at one level
				if(rxn in pathways_dict):
					print(pwy,subpwy,rxn)

				if(rxn not in pwy_data['reactions']):
					pwy_data['reactions'].append(rxn)
			
		mpfh.write("\t".join([pwy, pwy_data['name'],  "|".join(pwy_data['reactions']), "|".join(pwy_data['parent'])])+"\n")
