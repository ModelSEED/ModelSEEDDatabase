#!/usr/bin/env python
import os, sys, re, copy, time
from BiochemPy import Compounds

if(len(sys.argv)!=2 or os.path.isfile(sys.argv[1]) is False):
    print("Takes one argument, the path to and including the compound attributes file")
    sys.exit()

attributes_file=sys.argv[1]

#For logging purposes
log_fh = open("../../../Biochemistry/curation/all_updates.log","a+")
curation_source = attributes_file.split('/')[-2]
time_str = time.strftime('%Y-%m-%d', time.gmtime(time.time()))

#Load compounds
compounds_helper = Compounds()
compounds_dict = compounds_helper.loadCompounds()

Headers=list()
updated_cpds=list()
with open(attributes_file) as fh:
    for line in fh.readlines():
        line=line.strip()
        if(len(Headers)==0):
            Headers=line.split('\t')
            continue

        cpd_attributes=dict()
        array=line.split('\t',len(Headers))
        for i in range(len(Headers)):
            cpd_attributes[Headers[i].lower()]=array[i]

        #Check for compound
        if('id' not in cpd_attributes): 
            print("Warning: compounds' 'id' attribute not present")
            sys.exit()

        if(cpd_attributes['id'] not in compounds_dict):
            print("Warning: compound id ("+cpd_attributes['id']+") not in biochemistry database")
            continue

        #Iterating through attributes to assign them
        for attribute in cpd_attributes:
            if(attribute == 'id'):
                continue

            cpd_id = cpd_attributes['id']
            if(attribute not in compounds_dict[cpd_id]):
                print("Warning: compound attribute ("+attribute+") not recognized")
                continue

            if(cpd_attributes[attribute] != compounds_dict[cpd_id][attribute]):
                print("Updating " + cpd_id + ":" + attribute + \
                          " from " + str(compounds_dict[cpd_id][attribute]) + \
                          " to "+cpd_attributes[attribute])
                log_fh.write("\t".join([time_str, curation_source, cpd_id, attribute, \
                                            str(compounds_dict[cpd_id][attribute]), \
                                            cpd_attributes[attribute]])+"\n")

                compounds_dict[cpd_id][attribute] = cpd_attributes[attribute]

                if(cpd_id not in updated_cpds):
                    updated_cpds.append(cpd_id)

#        if(formula_charge_dict['formula'] != "null"):
#            compounds_dict[cpd]['formula']=formula_charge_dict['formula']
#            compounds_dict[cpd]['charge']=int(formula_charge_dict['charge'])

print("Saving "+str(len(updated_cpds))+" updated compounds")
compounds_helper.saveCompounds(compounds_dict)
log_fh.close()
