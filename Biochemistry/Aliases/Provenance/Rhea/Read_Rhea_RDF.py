#!/usr/bin/env python
import rdflib
from rdflib import Graph

graph = Graph()
graph.parse('Data/rhea.rdf')

coefficients_dict=dict()
reaction_participants=dict()
reactive_parts=dict()
participant_ids=dict()
chebi_ids=dict()
generic_ids=dict()
polymer_ids=dict()
compound_names=dict()
child_reactions=dict()
for subject, predicate, object in graph.triples((None, None, None)):
    if(predicate.endswith('/coefficient')):
        coefficient=str(object)
        #Fix issues with polymers
        if(coefficient == "2n" or coefficient == "Nplus1"):
            coefficient = "2"
        elif(coefficient == "Nminus1"):
            coefficient = "0"
        elif(coefficient == "N"):
            coefficient = "1"

        coefficients_dict[str(subject)]=coefficient

    if(predicate.endswith('directionalReaction')):
        parent_reaction = str(subject).split('/')[-1]
        child_reaction = str(object).split('/')[-1]
        if(parent_reaction not in child_reactions):
            child_reactions[parent_reaction]=list()
        child_reactions[parent_reaction].append(child_reaction)
        
    if((subject.endswith("_L") or subject.endswith("_R")) and "contains" in predicate):
        (reaction,side) = subject.rsplit('_')
        if(reaction not in reaction_participants):
            reaction_participants[reaction]=dict()
        if(side not in reaction_participants[reaction]):
            reaction_participants[reaction][side]=list()

        #here object is participant and predicate is coefficient
        reaction_participants[reaction][side].append([str(object),str(predicate)])

    if('accession' in predicate and 'GENERIC' in object):
        generic_ids[str(subject)]=object.split('_')[-1]
    
    if('accession' in predicate and 'POLYMER' in object):
        polymer_ids[str(subject)]=object.split('_')[-1]

    if('Compound' in subject and 'name' in predicate):
        if(str(subject) not in compound_names):
            compound_names[str(subject)]=list()
        compound_names[str(subject)].append(str(object))
        
    if("Compound" in object):
        participant_ids[str(subject)]=str(object)

    if("reactivePart" in predicate):
        reactive_parts[str(object)]=str(subject)

    #print("------------------------------")
    #print("Subject: ",subject)
    #print("Predicate: ",predicate)
    #print("Object: ",object)

    if("chebi" in predicate):
        chebi_ids[str(subject)]=object.split('_')[-1]

with open('ChEBI_rdf_names.tsv','w') as cpd_fh:
    cpd_fh.write('ID\tNAMES\n')
    for cpd in sorted(chebi_ids.keys()):
        cpd_fh.write(chebi_ids[cpd]+"\t"+"|".join(compound_names[cpd])+"\n")

with open('Rhea_rdf_names.tsv','w') as gp_fh:
    gp_fh.write('ID\tNAMES\n')
    for cpd in sorted(generic_ids.keys()):
        gp_fh.write(generic_ids[cpd]+"\t"+"|".join(compound_names[cpd])+"\n")
    for cpd in sorted(polymer_ids.keys()):
        gp_fh.write(polymer_ids[cpd]+"\t"+"|".join(compound_names[cpd])+"\n")

ofh = open("Rhea_reactions.tsv",'w')
efh = open("missing_rhea_entities.tsv",'w')
ofh.write("ID\tEQUATION\n")
for reaction in sorted(reaction_participants):
    reaction_list=list()
    complete_reaction=True
    for side in ('L','R'):
        side_list=list()
        for participant in reaction_participants[reaction][side]:
            (coefficient,compound,index)=(None,None,"0")
            if(participant[1] in coefficients_dict):
                coefficient = coefficients_dict[participant[1]]
            if(participant[0] in participant_ids):
                compound = ""
                if(participant_ids[participant[0]] in chebi_ids):
                    compound = chebi_ids[participant_ids[participant[0]]]
                elif(participant_ids[participant[0]] in generic_ids):
                    compound = generic_ids[participant_ids[participant[0]]]
                elif(participant_ids[participant[0]] in polymer_ids):
                    compound = polymer_ids[participant_ids[participant[0]]]
                else:
                    efh.write("Warning One: "+participant[0]+"\t"+participant_ids[participant[0]]+"\n")
                    complete_reaction = False
            else:
                efh.write("Warning Two: "+participant[0]+"\n")
                complete_reaction = False
            if(participant[0].endswith("_out")):
                index="1"
            if(coefficient is not None and compound is not None):
                side_list.append("("+coefficient+") "+compound+"["+index+"]")
         
        side_string=" + ".join(side_list)
        reaction_list.append(side_string)

    reaction_id=reaction.split('/')[-1]
    if(complete_reaction is True):
        ofh.write("\t".join([reaction_id," <=> ".join(reaction_list)])+"\n")
        if(reaction_id in child_reactions):
            for child_reaction in sorted(child_reactions[reaction_id]):
                ofh.write("\t".join([child_reaction," <=> ".join(reaction_list)])+"\n")
    else:
        efh.write("\t".join([reaction_id," <=> ".join(reaction_list)])+"\n")
