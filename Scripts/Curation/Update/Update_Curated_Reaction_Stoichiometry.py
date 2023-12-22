#!/usr/bin/env python
import os, sys, re, copy, time
from BiochemPy import Reactions

if(len(sys.argv)!=2 or os.path.isfile(sys.argv[1]) is False):
    print("Takes one argument, the path to and including the reaction stoichiometry file")
    sys.exit()

stoich_file=sys.argv[1]

#For logging purposes
log_fh = open("../../../Biochemistry/curation/all_updates.log","a+")
curation_source = stoich_file.split('/')[-2]
time_str = time.strftime('%Y-%m-%d', time.gmtime(time.time()))

#Load reactions
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

Headers=list()
updated_rxns=list()
with open(stoich_file) as fh:
    for line in fh.readlines():
        line=line.strip()

        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        #Enforced headers: 0:Reaction, 1:Compound, 2:Stoichiometry, 3:Compartment
        (line_reaction,line_compound,line_stoich,line_compartment)=line.split('\t')
        if(line_reaction not in reactions_dict):
            print("Skipping "+line_reaction+" as identifier not found in database")
            continue

        rxn_cpds_array = reactions_helper.parseStoich(reactions_dict[line_reaction]['stoichiometry'])
        updated_rxn=False
        for cpd_dict in rxn_cpds_array:
            if(cpd_dict['compound'] == line_compound):
                
#                print(line_reaction,cpd_dict['compound'],cpd_dict['coefficient'],cpd_dict['compartment'])

                if(cpd_dict['coefficient'] != int(line_stoich)):
                    print("Updating " + line_reaction + ":" + line_compound + " coefficient " + \
                              " from " + str(cpd_dict['coefficient']) + \
                              " to "+line_stoich)

                    log_fh.write("\t".join([time_str, curation_source, \
                                                line_reaction, line_compound, "coefficient", \
                                                str(cpd_dict['coefficient']), \
                                                line_stoich+"\n"]))

                    cpd_dict['coefficient'] = int(line_stoich)
                    updated_rxn=True

                if(cpd_dict['compartment'] != int(line_compartment)):
                    print("Updating " + line_reaction + ":" + line_compound + " compartment " + \
                              " from " + str(cpd_dict['compartment']) + \
                              " to "+line_compartment)

                    log_fh.write("\t".join([time_str, curation_source, \
                                                line_reaction, line_compound, "compartment", \
                                                str(cpd_dict['compartment']), \
                                                line_compartment+"\n"]))

                    cpd_dict['compartment'] = int(line_compartment)
                    updated_rxn=True
                
        if(updated_rxn is True and line_reaction not in updated_rxns):
            updated_rxns.append(line_reaction)

        stoichiometry = reactions_helper.buildStoich(rxn_cpds_array)
        reactions_dict[line_reaction]['stoichiometry'] = stoichiometry

print("Saving "+str(len(updated_rxns))+" updated reactions")
reactions_helper.saveReactions(reactions_dict)
log_fh.close()
