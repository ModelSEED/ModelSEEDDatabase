#!/usr/bin/env python
import os, sys, re, copy, time
from BiochemPy import Reactions

if(len(sys.argv)!=2 or os.path.isfile(sys.argv[1]) is False):
    print("Takes one argument, the path to and including the reaction stoichiometry file")
    sys.exit()

stoich_file=sys.argv[1]

#For logging purposes
log_fh = open("../../../Biochemistry/Curation/all_updates.log","a+")
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

        #Enforced headers: 0:Reaction, 1:Old Compound, 2:New Compound
        print(line)
        (line_reaction,line_old_compound,line_new_compound)=line.split('\t')
        if(line_reaction not in reactions_dict):
            print("Skipping "+line_reaction+" as identifier not found in database")
            continue

        updated_rxn=False
        rxn_cpds_array = reactions_dict[line_reaction]["stoichiometry"]
        stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)
            
        reactions_helper.replaceCompound(rxn_cpds_array,line_old_compound,line_new_compound)
        new_stoichiometry=reactions_helper.buildStoich(rxn_cpds_array)

        if(stoichiometry != new_stoichiometry):
            print("Replacting "+line_old_compound+" with "+line_new_compound+" in "+line_reaction)
            log_fh.write("\t".join([time_str, curation_source, \
                                    line_reaction, line_old_compound, line_new_compound+"\n"]))

            updated_rxn=True
                
        if(updated_rxn is True and line_reaction not in updated_rxns):
            updated_rxns.append(line_reaction)

if(len(updated_rxns)>0):
    print("Saving "+str(len(updated_rxns))+" updated reactions")
    reactions_helper.saveReactions(reactions_dict)
log_fh.close()
