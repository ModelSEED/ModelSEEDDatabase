#!/usr/bin/env python
import re, sys
#sys.path.append('../../../Libs/Python/')
#from BiochemPy import Compounds

pathways_dict = dict()

# Need to capture pathway hierarchy
path_parent_dict=dict()

###########################################################
# Parse modules first
# The pathway file has a lot of organism-specific
# pathways with no reactions and missing module links
# I want to mine the module file first to make it
# straightforward to skip the pathways
###########################################################

module=None
module_file = "data/module"
with open(module_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_lst=re.split('\s+',line)
		data=" ".join(tmp_lst[1:])
		if(tmp_lst[0] != ""):
			field = tmp_lst[0]

		if(field == "ENTRY"):
			module = tmp_lst[1]
			pathways_dict[module] = {'name':"",'description':"",
						 'reactions':[],
						 'parents':[]}

		elif(field == "NAME"):
			data=data.rstrip(';')
			pathways_dict[module]['name']=data
		elif(field == "DESCRIPTION"):
			pathways_dict[module]['description']=data
		elif(field == "REACTION"):
			data = data.split(' ')[0]
			rxns=re.split('[\s,]',data)
			for rxn in rxns:
				pathways_dict[module]['reactions'].append(rxn)
		elif(field == "PATHWAY"):
			pwy=data.split(" ")[0]
			pwy = pwy.replace('map','rn')
			if(pwy not in pathways_dict[module]['parents']):
				pathways_dict[module]['parents'].append(pwy)
		elif(field == "CLASS"):
			classes = data.split("; ")
			pathways_dict[module]['classes']=classes

		else:
			#print(field)
			pass

		if(field == "///"):

			if(len(pathways_dict[module]['reactions']) > 0):
				classes = pathways_dict[module]['classes']
				cls=classes[-1]

				if(cls not in pathways_dict[module]['parents']):
					pathways_dict[module]['parents'].append(cls)

				for i in range(len(classes)-1,-1,-1):
					if(classes[i] == "Pathway modules"):
						continue

					if(classes[i] not in path_parent_dict):
						path_parent_dict[classes[i]]=list()

					if(i>0):
						if(classes[i-1] == "Pathway modules"):
							continue

						if(classes[i-1] not in path_parent_dict[classes[i]]):
							path_parent_dict[classes[i]].append(classes[i-1])
					else:
						if(None not in path_parent_dict[classes[i]]):
							path_parent_dict[classes[i]].append(None)

###########################################################
# Parse pathways
###########################################################

pathway_file = "data/pathway"
pathway_id=None

with open(pathway_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_lst=re.split('\s+',line)
		data=" ".join(tmp_lst[1:])
		if(tmp_lst[0] != ""):
			field = tmp_lst[0]

		if(field == "ENTRY"):
			pathway_id = tmp_lst[1]
			pathways_dict[pathway_id] = {'name':"",'description':"",
						  'reactions':[],'classes':"",
						  'parents':[],'children':[]}

		elif(field == "NAME"):
			data=data.rstrip(';')
			pathways_dict[pathway_id]['name']=data
		elif(field == "DESCRIPTION"):
			pathways_dict[pathway_id]['description']=data
		elif(field == "REACTION"):
			rxn=data.split(" ")[0]
			pathways_dict[pathway_id]['reactions'].append(rxn)
		elif(field == "MODULE"):
			mod=data.split(" ")[0]
			if(mod in pathways_dict):
				pathways_dict[pathway_id]['children'].append(mod)
		elif(field == "CLASS"):
			classes = data.split("; ")
			pathways_dict[pathway_id]['classes']=classes
		else:
			#print(field)
			pass

		if(field == "///"):
			
			if(len(pathways_dict[pathway_id]['reactions']) > 0):
				classes = pathways_dict[pathway_id]['classes']
				cls=classes[-1]

				if(cls not in pathways_dict[pathway_id]['parents']):
					pathways_dict[pathway_id]['parents'].append(cls)

				for i in range(len(classes)-1,-1,-1):
					if(classes[i] not in path_parent_dict):
						path_parent_dict[classes[i]]=list()

					if(i>0):
						if(classes[i-1] not in path_parent_dict[classes[i]]):
							path_parent_dict[classes[i]].append(classes[i-1])
					else:
						if(None not in path_parent_dict[classes[i]]):
							path_parent_dict[classes[i]].append(None)

pwys_to_del = list()
for pwy in pathways_dict:
	if(len(pathways_dict[pwy]['reactions']) == 0):
		# don't remove pathways that have a list of modules
		if('children' in pathways_dict[pwy] and len(pathways_dict[pwy]['children'])>0):
			continue

		# most modules have reactions some don't but are still metabolic?
		if(pwy.startswith('M')):
			if("Pathway modules" in pathways_dict[pwy]['classes']):
				continue

			if("Metabolic capacity" in pathways_dict[pwy]['classes']):
				continue

		pwys_to_del.append(pwy)

for pwy in pwys_to_del:
	del(pathways_dict[pwy])

for pwy in pathways_dict:
	for parent in pathways_dict[pwy]['parents']:
		if(parent in pwys_to_del):
			pathways_dict[pwy]['parents'].remove(parent)

for child in path_parent_dict:
	if(child not in pathways_dict):
		pathways_dict[child]={'name':"",
				      'description':"",
				      'reactions':[],
				      'pathways':[],
				      'parents':[]}

	for parent in path_parent_dict[child]:
		if(parent==child):
			continue

		if(parent is None):
			continue

		if(parent not in pathways_dict[child]['parents']):
			pathways_dict[child]['parents'].append(parent)

		if(parent not in pathways_dict):
			pathways_dict[parent]={'name':"",
					       'description':"",
					       'reactions':[],
					       'pathways':[],
					       'parents':[]}

with open("KEGG_pathways.tsv",'w') as kpfh:
	kpfh.write("\t".join(["id","name","description","reactions","parent"])+"\n")
	for pwy in pathways_dict:
		kpfh.write("\t".join([pwy, pathways_dict[pwy]['name'], pathways_dict[pwy]['description'], "|".join(pathways_dict[pwy]['reactions']), "|".join(pathways_dict[pwy]['parents'])])+"\n")
