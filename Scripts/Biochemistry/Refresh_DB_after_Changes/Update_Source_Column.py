#!/usr/bin/env python
import sys
import os

sys.path.append('../../../Libs/Python/')
from BiochemPy import Reactions, Compounds

biochem_root = os.path.dirname(__file__)+'/../../../Biochemistry/'
curation_path = os.path.join(biochem_root,"Curation")
curators_list=list()
for dir in os.listdir(curation_path):
	if(os.path.isdir(biochem_root+'Curation/'+dir)):
		curators_list.append(dir)

source_classes_dict = dict()
with open('../../../Biochemistry/Aliases/Source_Classifiers.txt') as sc_fh:
	for line in sc_fh.readlines():
		line=line.strip('\r\n')
		(source,classifier) = line.split('\t')
		if(source in source_classes_dict):
			print("Warning: source listed more than once")
		source_classes_dict[source]=classifier

compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()
aliases_dict = compounds_helper.loadMSAliases()

updated_compounds_list=list()
for cpd in compounds_dict.keys():
	
	main_source_type = 'Orphan'

	if(cpd not in aliases_dict):
		#Set as Orphan
		if(compounds_dict[cpd]['source']!='Orphan'):
			compounds_dict[cpd]['source']='Orphan'
			print("New Orphan Compound:",cpd)
		continue

	for source_type in 'Primary Database', 'Secondary Database', 'Published Model', 'User':
		for source in aliases_dict[cpd]:
			if(source not in source_classes_dict and source in curators_list):
				main_source_type = "User"
				break
			if(source_classes_dict[source] == source_type):
				main_source_type = source_type
				break
		if(main_source_type != 'Orphan'):
				break
		
	if(compounds_dict[cpd]['source'].strip() == ''):
		print("New Source for",cpd,":",main_source_type)
		compounds_dict[cpd]['source']=main_source_type
		updated_compounds_list.append(cpd)
	elif(main_source_type != compounds_dict[cpd]['source']):
		print("Change in Source for "+cpd+": From",compounds_dict[cpd]['source'],"to",main_source_type)
		compounds_dict[cpd]['source']=main_source_type
		updated_compounds_list.append(cpd)

if(len(updated_compounds_list)>0):
	print("Updating source for "+str(len(updated_compounds_list))+" compounds")
	compounds_helper.saveCompounds(compounds_dict)

reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()
aliases_dict = reactions_helper.loadMSAliases()

updated_reactions_list=list()
for rxn in reactions_dict.keys():
	
	main_source_type = 'Orphan'

	if(rxn not in aliases_dict):
		#Set as Orphan
		if(reactions_dict[rxn]['source']!='Orphan'):
			reactions_dict[rxn]['source']='Orphan'
			print("New Orphan reaction:",rxn)
		continue

	for source_type in 'Primary Database', 'Secondary Database', 'Published Model', 'User':
		for source in aliases_dict[rxn]:
			if(source == 'rhea'):
				continue
			if(source not in source_classes_dict and source in curators_list):
				main_source_type = "User"
				break
			if(source_classes_dict[source] == source_type):
				main_source_type = source_type
				break
		if(main_source_type != 'Orphan'):
				break
		
	if(reactions_dict[rxn]['source'].strip() == ''):
		print("New Source for",rxn,":",main_source_type)
		reactions_dict[rxn]['source']=main_source_type
		updated_reactions_list.append(rxn)
	elif(main_source_type != reactions_dict[rxn]['source']):
		print("Change in Source for "+rxn+": From",reactions_dict[rxn]['source'],"to",main_source_type)
		reactions_dict[rxn]['source']=main_source_type
		updated_reactions_list.append(rxn)

if(len(updated_reactions_list)>0):
	print("Updating source for "+str(len(updated_reactions_list))+" reactions")
	reactions_helper.saveReactions(reactions_dict)
