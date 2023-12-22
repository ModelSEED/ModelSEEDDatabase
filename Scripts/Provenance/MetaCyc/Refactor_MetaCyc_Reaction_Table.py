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

###########################################################
# Load compounds to check
###########################################################

cpds_dict=dict()
with open('MetaCyc_compounds.tsv') as fh:
	for line in fh.readlines():
		line=line.strip('\r\n')
		cpd_id=line.split('\t')[0]
		cpds_dict[cpd_id]=line

class_dict = dict()
classes_dict = dict()
with open('data/classes.dat') as cfh:
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

			class_dict = {'id':data,'names':[],
							'formula':"null",
							'mass':"10000000",
							'charge':"10000000",
							'inchikey':'','smiles':''}

		elif(field == "SYNONYMS" or field.endswith('NAME')):
			if(field == "STRAIN-NAME" or field == "N-NAME" or field == "N+1-NAME" or field == "N-1-NAME"):
				continue

			name = convert_html_name(data)
			if(name not in class_dict['names']):
				class_dict['names'].append(name)
		elif(field.startswith("//")):
			class_line = "\t".join([class_dict['id'],"|".join(class_dict['names']),class_dict['formula'], \
								class_dict['mass'],class_dict['charge'],class_dict['inchikey'],class_dict['smiles']])
			classes_dict[class_dict['id']]=class_line

# There is a special compound that represents the electron and is missing
# From the compounds and classes file
electron_dict = {'id':'E-','names':['electron'],
					'mass':"0",'charge':"-1",
					'formula':"null",
					'inchikey':'','smiles':''}
electron_line = "\t".join([electron_dict['id'],"|".join(electron_dict['names']),electron_dict['formula'], \
					electron_dict['mass'],electron_dict['charge'],electron_dict['inchikey'],electron_dict['smiles']])
classes_dict['E-']=electron_line

reactions_list = list()
reaction_dict = dict()
rqd_classes_list = ['E-']
reagents = list()
with open('data/reactions.dat') as rfh:
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

			reaction_dict = {'id':data,
							'names':[],'ecs':[],
							'direction':None,
							'equation':""}

			# additional means of determining transport
			rxn_location = None
			reagents = list()

		elif(field == "SYNONYMS" or field.endswith('NAME')):

			name = convert_html_name(data)
			if(name not in reaction_dict['names']):
				reaction_dict['names'].append(name)

		elif(field == "DBLINKS"):
			data=re.sub('[\(\)]','',data)

			tmp_lst = data.split(' ')
			tmp_lst[1]=re.sub('"','',tmp_lst[1])
			if(tmp_lst[0] == "LIGAND-RXN"):
				print(tmp_lst[1])
			if(tmp_lst[0] == "RHEA"):
				print(tmp_lst[1])
		elif(field == "LEFT"):

			rgt = {'id':data,'cpt':'0','cft':'1','side':'L'}
			reagents.append(rgt)

		elif(field == "RIGHT"):

			rgt = {'id':data,'cpt':'0','cft':'1','side':'R'}
			reagents.append(rgt)

		elif(field == "^COEFFICIENT"):

			#Handling weird variations in coefficients
			data=re.sub(' ','',data)
			data=re.sub('\(','',data)
			data=re.sub('\)','',data)
			data=data.lower()
			coefficient = convert_html_name(data)
			reagents[-1]['cft']=coefficient

		elif(field == "^COMPARTMENT"):

			reagents[-1]['cpt']=data

		elif(field == "EC-NUMBER"):

			data=re.sub('\|','',data)
			ec_number=re.sub('^EC-','',data)
			if(len(ec_number.split('.')) == 2):
				ec_number+='.-.-'
			if(len(ec_number.split('.')) == 3):
				ec_number+='.-'
			if(ec_number not in reaction_dict['ecs']):
				reaction_dict['ecs'].append(ec_number)

		elif(field == "REACTION-DIRECTION"):

			if('LEFT-TO-RIGHT' in data):
				reaction_dict['direction']='=>'
			if('RIGHT-TO-LEFT' in data):
				reaction_dict['direction']='<='
			if('REVERSIBLE' in data):
				reaction_dict['direction']='<=>'

		elif(field == "RXN-LOCATIONS"):

			#Transport can be determined by reaction locations with two locations
			#But may need to determine which way from reagents?
			rxn_location=data

		elif(field.startswith('//')):

			rgts=list()
			pdts=list()
			update_cft = False
			update_cpt_list = list()
			reaction_is_complete=True
			for rgt in reagents:
				# Need to check that reagent is in compound list
				if(rgt['id'] not in cpds_dict):
					if(rgt['id'] in classes_dict and rgt['id'] not in rqd_classes_list):
						rqd_classes_list.append(rgt['id'])
					if(rgt['id'] not in classes_dict):
						print("ERROR: ",reaction_dict['id']+" contains unknown compound "+rgt['id'])
						reaction_is_complete=False

				# Need to check if the coefficients or compartments
				# should be standardized
				if(re.search('\D',rgt['cft']) is not None):
					update_cft = True

				if(rgt['cpt'] not in update_cpt_list):
					update_cpt_list.append(rgt['cpt'])

				if(rgt['side']=='L'):
					rgts.append(rgt)
				if(rgt['side']=='R'):
					pdts.append(rgt)

			if( len(rgts)==0 or len(pdts)==0 ):
				#print("Warning: Reaction",reaction_dict['id'],"doesn't have enough reagents")
				continue

			# fix missing direction
			if(reaction_dict['direction'] is None):
				reaction_dict['direction'] = '<=>'

			########################################
			# build equation
			########################################
			eqn = list()
			rgts_str = list()
			for rgt in rgts:
				rgt_str = "("+rgt['cft']+") "+rgt['id']+"["+rgt['cpt']+"]"
				rgts_str.append(rgt_str)
			eqn.append(" + ".join(rgts_str))
			eqn.append(reaction_dict['direction'])
			pdts_str = list()
			for pdt in pdts:
				pdt_str = "("+pdt['cft']+") "+pdt['id']+"["+pdt['cpt']+"]"
				pdts_str.append(pdt_str)
			eqn.append(" + ".join(pdts_str))
			eqn_str = " ".join(eqn)

			########################################
			# update compartments
			########################################
			
			# if they're all the same, there's no transport
			# but we do the replacement here to force them
			# to always be zero
			if(len(update_cpt_list) == 1):
				
				eqn_str = re.sub(update_cpt_list[0],'0',eqn_str)

				if(rxn_location is not None and len(re.findall('CCO',rxn_location))>1):
					#print("Warning: Reaction",reaction_dict['id'],"should have more than one compartment")
					pass

			elif(len(update_cpt_list)==2 and '0' in update_cpt_list):

				# transport reaction that has either CCO-IN or CCO-OUT
				# we declare CCO-IN to be cpt 'zero' and CCO-OUT to be cpt 'one'
				if('CCO-IN' in update_cpt_list):
					# replace [0] with [1]
					eqn_str = re.sub('\[0\]','[1]',eqn_str)
					# then replace 'CCO-IN' with 0
					eqn_str = re.sub('CCO-IN','0',eqn_str)
				elif('CCO-OUT' in update_cpt_list):
					# replace 'CCO-OUT' with 1
					eqn_str = re.sub('CCO-OUT','1',eqn_str)

			elif(len(update_cpt_list)>=2):
				
				# start with CCO-IN as 0
				if('CCO-IN' in update_cpt_list):
					# if CCO-IN is 0 then the other cpt is 1
					for cpt in update_cpt_list:
						if(cpt == 'CCO-IN'):
							eqn_str = re.sub(cpt,'0',eqn_str)
						# special rule #1 if both CCO-IN and CCO-CYTOSOL present
						# in a multi-compartment equation then
						# both are assigned 0
						elif(cpt == 'CCO-CYTOSOL' and len(update_cpt_list)>2):
							eqn_str = re.sub(cpt,'0',eqn_str)
						# any other compartment except 0 is then assigned as 1
						# in the case of multiple compartments, including
						# CCO-MEMBRANE and CCO-MIDDLE, these are also assigned 1
						elif(cpt != '0'):
							eqn_str = re.sub(cpt,'1',eqn_str)

				else:
					# as a reminder, many equations contained 'CCO-IN'
					# then have already been processed. There are a few
					# reactions that have 'CCO-OUT' and others instead

					# if CCO-OUT is 1 then the other cpt is 0	
					for cpt in update_cpt_list:
						if(cpt == 'CCO-OUT'):
							eqn_str = re.sub(cpt,'1',eqn_str)
						else:
							eqn_str = re.sub(cpt,'0',eqn_str)

			########################################
			# update coefficients
			########################################

			if(update_cft is True):

				# this is a fix for a single reaction (RXN-12218)
				# that appears to have been mis-annotated. It'll
				# end up unbalanced anyway
				eqn_str = re.sub('m\+q','1',eqn_str)

				# find common occurences
				minus_cft_tuples = re.findall('(([nm])-(\d))\)',eqn_str)
				plus_cft_tuples = re.findall('(([nm])\+(\d))\)',eqn_str)
				multi_cft_tuples = re.findall('((\d)([nm]))\)',eqn_str)

				# first we handle deductions, so far only one type
				# is done per equation, i.e. either n-1 or m-1
				if(len(minus_cft_tuples)>0):
					# we go through these and count how many
					# this is important as there may be multiple
					# instances in a single equation
					minus_cft_variable = 'n' # this may change
					minus_cft_count = 0
					for entry in minus_cft_tuples:
						if(entry[1] != minus_cft_variable):
							# i.e. if it's not 'n', sometimes it's 'm'
							minus_cft_variable = entry[1]
						minus_cft_count += int(entry[2])
						# after counting each one, we replace the instance
						# with the default value of 1
						eqn_str = re.sub(re.escape(entry[0])+'\)','1)',eqn_str)
					# finally, we replace the original stand alone variable
					# with the count. So, for example, this means that we
					# change: n -> n-1 + n-1 to 2 -> 1 + 1
					eqn_str = re.sub('\('+re.escape(minus_cft_variable)+'\)','('+str(minus_cft_count)+')',eqn_str)
				
				# re-running this because the previous clause may result in changing
				# m+n-1 to m+1 so we get new instances to fix
				plus_cft_tuples = re.findall('(([nm])\+(\d))\)',eqn_str)
				# secondly we handle additions, usually n+1 or m+1
				# but in some cases the integer is higher
				if(len(plus_cft_tuples)>0):
					# we go through these and count how many
					# this is important as there may be multiple
					# instances in a single equation
					plus_cft_variable = 'n' # this may change
					plus_cft_count = 0
					for entry in plus_cft_tuples:
						if(entry[1] != plus_cft_variable):
							# i.e. if it's not 'n', sometimes it's 'm'
							minus_cft_variable = entry[1]
						plus_cft_count += int(entry[2])
						# after counting each one, we replace the instance
						# with the retrieved value incremented by the default value of 1
						eqn_str = re.sub(re.escape(entry[0])+'\)',str(1+plus_cft_count)+')',eqn_str)
				
				# thirdly, we handle multiplication, usually 2n,3n etc.
				# as n defaults to 1, this is simple to fix
				if(len(multi_cft_tuples)>0):
					for entry in multi_cft_tuples:
						eqn_str = re.sub('\('+entry[0]+'\)','('+str(entry[1])+')',eqn_str)
						
				# finally, we replace any remaining stand alone variable
				# with the 1. So, for example, this means that we
				# change: n -> 1
				if(re.search('\(\D\)',eqn_str) is not None):
					for entry in re.findall('\(\D\)',eqn_str):
						eqn_str = re.sub(re.escape(entry),'(1)',eqn_str)

			reaction_dict['equation'] = eqn_str
			if(reaction_is_complete is True):
				reactions_list.append(reaction_dict)

with open("MetaCyc_reactions.tsv",'w') as mrfh:
	mrfh.write("\t".join(["id","names","equation","enzymes"])+"\n")
	for rxn in reactions_list:
		mrfh.write("\t".join([rxn['id'],"|".join(rxn['names']),rxn['equation'], "|".join(rxn['ecs'])])+"\n")

with open("MetaCyc_compounds.tsv",'a') as mcfh:
	for rqd in sorted(rqd_classes_list):
		mcfh.write(classes_dict[rqd]+'\n')
