#!/usr/bin/env python
import re, sys
#sys.path.append('../../../Libs/Python/')
#from BiochemPy import Compounds

def convert_KEGG_equation(equation):

	###########################################################
	# Account for all variables and their modifiers
	###########################################################
	#Two possible variables (n and m)
	for variable_group in ('n','m'):
		# Situation 1: minus symbol
		# Assume variable is 2 and replace 'variable-1' with 1
		if(variable_group+'-1' in equation):
			equation = re.sub(variable_group+'-1','1',equation)
			equation = re.sub(variable_group,'2',equation)

		# Situation 2: plus symbol (variable+)
		# Assume variable is 1 and replace 'variable+n' with 1+n
		if(variable_group+'+' in equation):
			match_list = re.findall(variable_group+'\+(\d)',equation)
			for match in match_list:
				match_str=variable_group+'\+'+match
				variable_mod=1+int(match)
				equation = re.sub(match_str,str(variable_mod),equation,count=1)
			
		# Situation 3: plus symbol (+variable)
		# Arises if (n+m) is encountered in Situation 2
		# Assume variable is 1 and replace ('1+variable') with 2
		if('1+'+variable_group in equation):
			equation = re.sub('1\+'+variable_group,'2',equation)
	
		# Situation 4: variable with multiplier (n*variable)
		# Assume variable is 1 and replace (n*variable) with n
		if(variable_group in equation):
			match_list = re.findall('(\d)'+variable_group,equation)
			for match in match_list:
				match_str=match+variable_group
				variable_mod=int(match)*1
				equation = re.sub(match_str,str(variable_mod),equation,count=1)
		
		# Situation 5: all remaining instances of variable
		# Assume variable is 1
		if(variable_group in equation):
			equation = re.sub(variable_group,'1',equation)

	###########################################################
	# Reorganize equation
	###########################################################

	# Split equation
	rgts_list=equation.split('<=>')
	for h in range(len(rgts_list)):
		rgt_str = rgts_list[h]

		# Split reagents
		rgt_list = rgt_str.split('+')
		for i in range(len(rgt_list)):

			# Strip whitespace
			rgt_list[i]=rgt_list[i].strip()
		
			# Normalize position of bracketed stoichiometries
			match_mod = re.search('\(\d\)$',rgt_list[i])
			if(match_mod is not None):
				match_str = match_mod.group(0)
				rgt_list[i] = re.sub(re.escape(match_str),'',rgt_list[i])
				rgt_list[i] = match_str+' '+rgt_list[i]

			# Detecting stoichiometry > 1 for each reagent
			stoich_cpd_list=rgt_list[i].split(' ')
			if(len(stoich_cpd_list)>1):

				# Ignore (by force-removing) all stoichiometries of 1
				if(stoich_cpd_list[0] == '1' or stoich_cpd_list[0] == '(1)'):
					rgt_list[i]=stoich_cpd_list[1]
				# Bracket stoichiometries correctly if not so
				elif('(' not in stoich_cpd_list[0] and ')' not in stoich_cpd_list[0]):
					rgt_list[i]='('+stoich_cpd_list[0]+') '+stoich_cpd_list[1]
				# Stoichiometries are already correctly bracketed
				else:
					#print(stoich_cpd_list[0])
					pass

			# Add compartment index (default to 0)
			rgt_list[i]+='[0]'

		# Build new reagents list
		new_rgts_str = ' + '.join(rgt_list)
		rgts_list[h]=new_rgts_str
	
	# Build new equation
	new_equation = ' <=> '.join(rgts_list)
	return new_equation


###########################################################
# Load compounds and glycans to check
###########################################################

cpds_dict=dict()
with open('KEGG_compounds.tsv') as fh:
	for line in fh.readlines():
		line=line.strip('\r\n')
		cpd_id=line.split('\t')[0]
		cpds_dict[cpd_id]=line

glcs_dict=dict()
with open('KEGG_glycans.tsv') as fh:
	for line in fh.readlines():
		line=line.strip('\r\n')
		glc_id=line.split('\t')[0]
		glcs_dict[glc_id]=line

# There is one missing glycan, and I'm manually entering it here
missing_glc_name='alpha-Kdo-(2->8)-[alpha-Kdo-(2->4)]-alpha-Kdo-(2->4)-alpha-Kdo-(2->6)-lipid IVA'
glcs_dict['G13061']="\t".join(['G13061',missing_glc_name,'null','10000000','10000000','',''])

###########################################################
# Parse reactions
###########################################################

reaction_file = "data/reaction"
reaction_dict = dict()
reactions_dict = list()
required_glcs_list = list()
with open(reaction_file) as cfh:
	for line in cfh.readlines():
		line=line.strip('\r\n')
		tmp_list=re.split('\s+',line)
		data=" ".join(tmp_list[1:])
		if(tmp_list[0] != ""):
			field = tmp_list[0]

		if(field == "ENTRY"):
			reaction_dict = {'id':tmp_list[1],
							'names':[],'ecs':[],
							'equation':"",
							'pathway':""}

		elif(field == "NAME"):
			data=data.rstrip(';')
			reaction_dict['names'].append(data)
		elif(field == "EQUATION"):
			equation = convert_KEGG_equation(data)
			reaction_dict['equation']=equation
		elif(field == "ENZYME"):
			for ec in data.split(' '):
				reaction_dict['ecs'].append(ec)
		else:
			#print(field)
			pass

		if(field == "///"):

			reaction_is_complete=True
			for entry in reaction_dict['equation'].split(' '):
				cpd=entry.replace('[0]','')
				if(cpd.startswith('C') and cpd not in cpds_dict):
					print("ERROR: ",reaction_dict['id']+" contains unknown compound "+cpd)
					reaction_is_complete=False
				if(cpd.startswith('G')):
					if(cpd not in glcs_dict):
						print("ERROR: ",reaction_dict['id']+" contains unknown glycan "+cpd)
						reaction_is_complete=False
					elif(cpd not in required_glcs_list):
						required_glcs_list.append(cpd)

			if(reaction_is_complete is True):
				reactions_dict.append(reaction_dict)

with open("KEGG_reactions.tsv",'w') as krfh:
	krfh.write("\t".join(["id","names","equation","enzymes"])+"\n")
	for rxn in reactions_dict:
		krfh.write("\t".join([rxn['id'],"|".join(rxn['names']),rxn['equation'], "|".join(rxn['ecs'])])+"\n")

with open("KEGG_compounds.tsv",'a') as kcfh:
	for glc in sorted(required_glcs_list):
		kcfh.write(glcs_dict[glc]+'\n')
