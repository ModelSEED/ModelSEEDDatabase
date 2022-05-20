import os
import re
import json
import copy
import glob
import itertools
from csv import DictReader

class Reactions:
    def __init__(self, biochem_root='../../../Biochemistry/',
                 rxns_file='reaction_00.tsv'):

        self.BiochemRoot = os.path.dirname(__file__)+'/'+biochem_root
        self.RxnsFile = self.BiochemRoot + rxns_file
        self.AliasFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
        self.NameFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Reaction_Names.txt"
        self.PwyFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Reaction_Pathways.txt"
        self.ECFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Reaction_ECs.txt"

        reader = DictReader(open(self.RxnsFile), dialect='excel-tab')
        self.Headers = reader.fieldnames

        from BiochemPy import Compounds
        self.CompoundsHelper = Compounds()
        self.Compounds_Dict = self.CompoundsHelper.loadCompounds()
            
    def loadReactions(self):

        search_path = os.path.join(self.BiochemRoot,"reaction_*.json")
        rxns_dict = dict()
        for reactions_file in sorted(glob.glob(search_path)):
            with open(reactions_file) as json_file_handle:
                rxns_list = json.load(json_file_handle)
                for rxn_obj in rxns_list:
                    for key in rxn_obj:
                        if(isinstance(rxn_obj[key],list)):
                            for i in range(len(rxn_obj[key])):
                                if(rxn_obj[key][i] is None):
                                    rxn_obj[key][i]="null"
                            
                        if(isinstance(rxn_obj[key],dict)):
                            for entry in rxn_obj[key]:
                                if(rxn_obj[key][entry] is None):
                                    rxn_obj[key][entry]="null"
                            
                        if(rxn_obj[key] is None):
                            rxn_obj[key]="null"
                        
                    rxns_dict[rxn_obj['id']]=rxn_obj
                    
        return rxns_dict

    def loadReactions_tsv(self):
        print("WARNING: This function is currently redundant and will only load one file!")
        reader = DictReader(open(self.RxnsFile), dialect='excel-tab')
        type_mapping = {"is_transport": int, "is_obsolete": int,
                        "deltag": float, "deltagerr": float}
        lists = ["aliases","pathways","ec_numbers","notes"]
        dicts = []

        rxns_dict = dict()
        for line in reader:
            for list_type in lists:
                if(line[list_type] != "null"):
                    line[list_type]=line[list_type].split("|")
            for dict_type in dicts:
                if(line[dict_type] != "null"):
                    entries = line[dict_type].split('|')
                    line[dict_type]=dict()
                    for entry in entries:
                        (type,list) = entry.split(':')
                        line[dict_type][type]=list
            for heading, target_type in type_mapping.items():
                try:
                    line[heading] = target_type(line[heading])
                except ValueError:  # Generally caused by "null" strings
                    line[heading] = None
            rxns_dict[line['id']] = line

        return rxns_dict

    def parseEquation(self, equation_string):
        rxn_cpds_array = list()
        reagent=-1
        coeff=1
        index=0
        for text in equation_string.split(" "):
            if(text == "+"):
                continue

            match=re.search('^<?[=-]+>?$',text)
            if(match is not None):
                reagent=1

            match=re.search('^\((\d+(?:\.\d+)?)\)$',text)
            if(match is not None):
                coeff=match.group(1)

                # Correct for redundant ".0" in floats
                coeff = float(coeff)
                if (str(coeff)[-2:] == ".0"):
                    coeff = int(round(coeff))

            match=re.search('^(cpd\d{5})\[(\d)\]$',text)
            if(match is not None):

                #Side of equation
                coeff=coeff*reagent

                (cpd,cpt)=(match.group(1),match.group(2))
                rgt_id = cpd + "_" + cpt + str(index)
                cpt = int(cpt)
                name = self.Compounds_Dict[cpd]["name"]
                formula = self.Compounds_Dict[cpd]["formula"]
                charge = self.Compounds_Dict[cpd]["charge"]

                rxn_cpds_array.append({"reagent": rgt_id, "coefficient": coeff,
                                       "compound": cpd, "compartment": cpt,
                                       "index": index, "name": name,
                                       "formula": formula, "charge": charge})

                #Need to reset coeff for next compound
                coeff=1

        return rxn_cpds_array

    def parseStoich(self, stoichiometry):
        rxn_cpds_array = list()

        #For empty reaction
        if(stoichiometry == ""):
            return rxn_cpds_array

        for rgt in stoichiometry.split(";"):
            (coeff, cpd, cpt, index, name) = rgt.split(":", 4)
            rgt_id = cpd + "_" + cpt + index

            coeff = float(coeff)

            # Correct for redundant ".0" in floats
            if (str(coeff)[-2:] == ".0"):
                coeff = int(round(coeff))

            cpt = int(cpt)
            index = int(index)

            rxn_cpds_array.append({"reagent": rgt_id, "coefficient": coeff,
                                   "compound": cpd, "compartment": cpt,
                                   "index": index, "name": name,
                                   "formula": self.Compounds_Dict[cpd][
                                       "formula"],
                                   "charge": self.Compounds_Dict[cpd][
                                       "charge"]})
        return rxn_cpds_array

    def parseStoichOnt(self, stoichiometry):
        rxn_cpds_dict = dict()

        #For empty reaction
        if(stoichiometry == ""):
            return rxn_cpds_array

        for rgt in stoichiometry.split(";"):
            (coeff, cpd, cpt, index, name) = rgt.split(":", 4)
            cpd_cpt_tuple = (cpd,cpt)
            rxn_cpds_dict[cpd_cpt_tuple]=coeff

        return rxn_cpds_dict

    # The basis for this code, and producing combinations of ontologically related reactions
    # was found in Filipe's code (see commit: 92db86)
    def generateOntologyReactionCodes(self, rxn_id, rxn_cpds, cpds_neighbors):

        # returns list of reaction codes to match with biochemistry
        new_codes = dict()

        replacements = list()
        for cpd_cpt_tuple in rxn_cpds:
            replace_list = list()
            cpd_id = cpd_cpt_tuple[0]
            if cpd_id in cpds_neighbors:
                for neighbor_id in cpds_neighbors[cpd_id]:
                    replace_list.append((cpd_id,neighbor_id))

            if len(replace_list) > 0:
                replacements.append(replace_list)

        # Iterate through different numbers of compounds to replace
        # i.e. replace 1 compound, replace 2 compounds etc.
        # The output is a list of all the possible combination of replacements to explore
        replacement_product=list()
        for n_cpds in range(1,len(replacements)+1):
            combination = list(itertools.combinations(replacements,n_cpds))
            for entry in combination:
                product_list = list(itertools.product(*entry))
                replacement_product += product_list

        if(len(replacements) == 0):
            return new_codes

        for entry in replacement_product:

            # Old code assumed that all "new" compounds were unique
            # cpd_swap_dict = {x:y for x, y in entry}
            # new_swapped_rxn_cpds = { (x if not x in cpd_swap_dict else cpd_swap_dict[x], c):y 
            #                          for (x, c), y in rxn_cpds.items() }

            # Regenerate array of cpd dicts for use with generateCode()
            swapped_rxn_cpds_array=list()
            for (cpd, cpt), coeff in rxn_cpds.items():
                new_cpd = cpd
                for old, new in entry:
                    if(cpd == old):
                        new_cpd = new
                reagent = { "reagent":new_cpd+'_'+cpt+'0',
                            "compartment":cpt,
                            "coefficient":float(coeff) }

                # Correct for redundant ".0" in floats
                if (str(reagent["coefficient"])[-2:] == ".0"):
                    reagent["coefficient"] = int(round(reagent["coefficient"]))

                swapped_rxn_cpds_array.append(reagent)
            new_code = self.generateCode(swapped_rxn_cpds_array)
            new_codes[new_code]=entry
        return new_codes

    @staticmethod
    def isTransport(rxn_cpds_array):
        compartments_dict=dict()
        for rgt in rxn_cpds_array:
            compartments_dict[rgt['compartment']]=1
        if(len(compartments_dict.keys())>1):
            return 1
        else:
            return 0

    def generateCodes(self, rxns_dict,check_obsolete=True):
        codes_dict=dict()
        for rxn in rxns_dict:
            if(rxns_dict[rxn]['status']=="EMPTY"):
                continue
            if(check_obsolete is False and rxns_dict[rxn]['is_obsolete']==1):
                continue
            rxn_cpds_array = self.parseStoich(rxns_dict[rxn]['stoichiometry'])
            code = self.generateCode(rxn_cpds_array)
            if(code not in codes_dict):
                codes_dict[code]=dict()
            codes_dict[code][rxn]=1
        return codes_dict

    def generateCode(self,rxn_cpds_array):

        #It matters if its a transport reaction, and we include protons when matching transport
        is_transport = self.isTransport(rxn_cpds_array)

        #It matters which side of the equation, so build reagents and products arrays
        reagents=list()
        products=list()
        for rgt in sorted(rxn_cpds_array, key=lambda x: ( x["reagent"], x["coefficient"] )):
            #skip protons
            if("cpd00067" in rgt["reagent"] and is_transport == 0):
                continue

            if(rgt["coefficient"]<0):
                reagents.append(rgt["reagent"]+":"+str(abs(rgt["coefficient"])))
            if(rgt["coefficient"]>0):
                products.append(rgt["reagent"]+":"+str(abs(rgt["coefficient"])))

        rgt_string = "|".join(reagents)
        pdt_string = "|".join(products)
        #Sorting the overall strings here helps with matching transporters
        rxn_string = "|=|".join(sorted([rgt_string,pdt_string]))
        return rxn_string

    @staticmethod
    def buildStoich(rxn_cpds_array):
        stoichiometry_array = list()
        for rgt in sorted(rxn_cpds_array, key=lambda x: (
                int(x["coefficient"] > 0), x["reagent"])):

            # Correct for redundant ".0" in floats
            if (str(rgt["coefficient"])[-2:] == ".0"):
                rgt["coefficient"] = int(round(rgt["coefficient"]))

            coeff = str(rgt["coefficient"])
            cpt = str(rgt["compartment"])
            idx = str(rgt["index"])

            rgt_string = ":".join(
                [coeff, rgt["compound"], cpt,
                 idx, '"'+rgt["name"]+'"'])
            stoichiometry_array.append(rgt_string)
        stoichiometry_string = ";".join(stoichiometry_array)
        return stoichiometry_string

    @staticmethod
    def removeCpdRedundancy(rgts_array):

        rgts_dict = dict()
        for rgt in rgts_array:
            if (rgt["reagent"] not in rgts_dict):
                rgts_dict[rgt["reagent"]] = 0
            rgts_dict[rgt["reagent"]] += float(rgt["coefficient"])

        new_rgts_array=list()
        for rgt in rgts_array:
            if (rgts_dict[rgt["reagent"]] == 0):
                continue

            rgt["coefficient"]=rgts_dict[rgt["reagent"]]

            # Correct for redundant ".0" in floats
            if (str(rgt["coefficient"])[-2:] == ".0"):
                rgt["coefficient"] = int(round(rgt["coefficient"]))

            new_rgts_array.append(rgt)
            
            #Trick to exclude reagent if it appears in array more than once
            rgts_dict[rgt["reagent"]]=0
            

        return new_rgts_array

    def balanceReaction(self, rgts_array):
        if (len(rgts_array) == 0):
            return "EMPTY"

        ########################################
        # Check that each reagent is either a 
        # different compound or in a different
        # compartment, and report.
        ########################################
        rgts_dict = dict()
        for rgt in rgts_array:
            if (rgt["reagent"] not in rgts_dict):
                rgts_dict[rgt["reagent"]] = 0
            rgts_dict[rgt["reagent"]] += 1

        for rgt in rgts_dict.keys():
            if (rgts_dict[rgt] > 1):
                return "Duplicate reagents"

        ########################################
        # Check for duplicate compounds in
        # different compartments, these are 
        # balanced directly.
        #######################################
        cpds_coeff_dict = dict()
        for rgt in rgts_array:
            cpd = rgt["compound"]
            if (cpd not in cpds_coeff_dict):
                cpds_coeff_dict[cpd] = 0

            # Use float() because you can get real coefficients
            cpds_coeff_dict[cpd] += float(rgt["coefficient"])

        # Build dict of compounds
        cpds_dict = dict()
        for rgt in rgts_array:
            #Skip trans-compartmental compounds
            if (cpds_coeff_dict[rgt["compound"]] == 0):
                continue

            proxy_rgt=copy.deepcopy(rgt)
            proxy_rgt["coefficient"] = cpds_coeff_dict[rgt["compound"]]
            cpds_dict[rgt["compound"]] = proxy_rgt

        ########################################
        # Check for duplicate elements, across
        # all compounds, these are balanced 
        # directly.
        #######################################
        rxn_net_charge = 0.0
        rxn_net_mass = dict()
        cpdformerror=list()
        for cpd in cpds_dict.keys():
            cpd_atoms = self.CompoundsHelper.parseFormula(
                cpds_dict[cpd]["formula"])

            if (len(cpd_atoms.keys()) == 0):
                #Here we can skip photons and electrons
                #They are the valid compounds with no mass
                if(cpd=='cpd11632' or cpd=='cpd12713'):
                    pass
                else:
                    cpdformerror.append(cpd)

            cpd_coeff_charge = float(cpds_dict[cpd]["charge"]) * float(
                cpds_dict[cpd]["coefficient"])
            rxn_net_charge += cpd_coeff_charge

            for atom in cpd_atoms.keys():
                atom_coeff_mass = float(cpd_atoms[atom]) * float(
                    cpds_dict[cpd]["coefficient"])

                if (atom not in rxn_net_mass.keys()):
                    rxn_net_mass[atom] = 0.0

                rxn_net_mass[atom] += atom_coeff_mass

        if(len(cpdformerror)>0):
            return "CPDFORMERROR"

        # Round out tiny numbers that occur because we add/substract floats
        # Threshold of 1e-6 found to capture all these instances without
        # removing actual small differences in mass.
        for atom in rxn_net_mass.keys():
            if (rxn_net_mass[atom] > -1e-6 and rxn_net_mass[atom] < 1e-6):
                rxn_net_mass[atom] = 0

        if (rxn_net_charge > -1e-6 and rxn_net_charge < 1e-6):
            rxn_net_charge = 0

        # Report any imbalance
        imbalanced_atoms_array = list()
        for atom in sorted(rxn_net_mass.keys()):
            if (rxn_net_mass[atom] == 0):
                continue

            rxn_net_mass[atom] = "{0:.2f}".format(rxn_net_mass[atom])

            # Correct for redundant ".00" in floats
            if (rxn_net_mass[atom][-3:] == ".00"):
                rxn_net_mass[atom] = str(int(float(rxn_net_mass[atom])))
    
            imbalanced_atoms_array.append(atom + ":" + rxn_net_mass[atom])

        rxn_net_charge = "{0:.2f}".format(rxn_net_charge)

        # Correct for redundant ".00" in floats
        if (rxn_net_charge[-3:] == ".00"):
            rxn_net_charge = str(int(float(rxn_net_charge)))

        status = ""

        if (len(imbalanced_atoms_array) > 0):
            status = "MI:" + "/".join(imbalanced_atoms_array)

        if (rxn_net_charge != "0"):
            if (len(status) == 0):
                status = "CI:" + rxn_net_charge
            else:
                status += "|CI:" + rxn_net_charge

        if (status == ""):
            status = "OK"
            
        return status

    def adjustCompound(self, rxn_cpds_array, compound, adjustment,
                       compartment=0):

        if (adjustment == 0):
            return rxn_cpds_array

        ######################################################################
        # We will always assume to adjust a compound automatically
        # in the compartment indexed as zero, unless otherwise specified.
        # This answers the question of how to handle transporters.
        ######################################################################

        # Check to see if it already exists
        cpd_exists = 0
        cpd_remove = {}
        for rgt in rxn_cpds_array:
            if (rgt["compound"] == compound and
                        rgt["compartment"] == compartment):
                rgt["coefficient"] -= adjustment
                cpd_exists = 1
                if(rgt["coefficient"] == 0):
                    cpd_remove=rgt

        if (cpd_exists != 1):
            rgt_id = compound + "_" + str(compartment) + "0"

            rxn_cpds_array.append(
                {"reagent": rgt_id, "coefficient": 0-adjustment,
                 "compound": compound, "compartment": compartment, "index": 0,
                 "name": self.Compounds_Dict[compound]["name"],
                 "formula": self.Compounds_Dict[compound]["formula"],
                 "charge": self.Compounds_Dict[compound]["charge"]})

        if(len(cpd_remove.keys())>0):
            rxn_cpds_array.remove(cpd_remove)

        #Got to adjust for floats
        for rgt in rxn_cpds_array:
            if (str(rgt["coefficient"])[-2:] == ".0"):
                rgt["coefficient"] = int(round(rgt["coefficient"]))

        return

    def replaceCompound(self, rxn_cpds_array, old_compound, new_compound):

        ######################################################################
        # We will always assume that we will maintain the coefficient.
        # We will always assume that we will replace in all compartments.
        # The adjustment will fail silently, returning an empty array
        # if the old_compound cannot be found.
        ######################################################################

        found_cpd=False
        for rgt in rxn_cpds_array:
            if (rgt["compound"] == old_compound):
                found_cpd=True
                rgt["compound"]=new_compound
                rgt["reagent"]=new_compound + "_" + str(rgt["compartment"]) + "0"
                rgt["name"]=self.Compounds_Dict[new_compound]['name']

        return found_cpd

    def rebuildReaction(self, reaction_dict, stoichiometry=None):
        # Retrieve/Assign stoich
        if(stoichiometry is None):
            stoichiometry = reaction_dict['stoichiometry']
        else:
            reaction_dict["stoichiometry"] = stoichiometry

        # Build list of "reagents" and "products"
        rxn_cpds_array = self.parseStoich(stoichiometry)
        reagents_array = list()
        products_array = list()
        compound_ids_dict = dict()
        for rgt in rxn_cpds_array:
            compound_ids_dict[rgt["compound"]] = 1
            if (rgt["coefficient"] > 0):
                products_array.append(rgt)
            else:
                reagents_array.append(rgt)

        rgts_str__array = list()
        for rgt in reagents_array:
            id_string = "(" + str(abs(rgt["coefficient"])) + ") " + rgt[
                "compound"] + "[" + str(rgt["compartment"]) + "]"
            rgts_str__array.append(id_string)

        equation_array = list()
        code_array = list()
        definition_array = list()

        equation_array.append(" + ".join(rgts_str__array))
        definition_array.append(" + ".join(rgts_str__array))
        code_array.append(
            " + ".join(x for x in rgts_str__array if "cpd00067" not in x))

        code_array.append("<=>")
        if (reaction_dict["direction"] == "="):
            equation_array.append("<=>")
            definition_array.append("<=>")
        elif (reaction_dict["direction"] == "<"):
            equation_array.append("<=")
            definition_array.append("<=")
        else:
            equation_array.append("=>")
            definition_array.append("=>")

        pdts_str_array = list()
        for rgt in products_array:
            id_string = "(" + str(abs(rgt["coefficient"])) + ") " + rgt[
                "compound"] + "[" + str(rgt["compartment"]) + "]"
            pdts_str_array.append(id_string)

        equation_array.append(" + ".join(pdts_str_array))
        definition_array.append(" + ".join(pdts_str_array))
        code_array.append(
            " + ".join(x for x in pdts_str_array if "cpd00067" not in x))

        reaction_dict["code"] = " ".join(code_array)
        reaction_dict["equation"] = " ".join(equation_array)
        reaction_dict["definition"] = " ".join(definition_array)
        reaction_dict["compound_ids"] = ";".join(
            sorted(compound_ids_dict.keys()))

        # Replace ids with names in Definition
        for cpd_id in compound_ids_dict.keys():
            if (cpd_id in reaction_dict["definition"]):
                reaction_dict["definition"] = reaction_dict[
                    "definition"].replace(cpd_id,
                                          self.Compounds_Dict[cpd_id]["name"])

        # Define if transport?

        return

    def saveECs(self, ecs_dict):
        ecs_root = os.path.splitext(self.ECFile)[0]

        # Print to TXT
        ecs_file = open(ecs_root + ".txt", 'w')
        ecs_file.write("\t".join(("ModelSEED ID","External ID","Source")) + "\n")
        for rxn in sorted(ecs_dict.keys()):
            for name in sorted(ecs_dict[rxn]):
                ecs_file.write("\t".join((rxn,name,'Enzyme Class')) + "\n")
        ecs_file.close()

    def saveNames(self, names_dict):
        names_root = os.path.splitext(self.NameFile)[0]

        # Print to TXT
        names_file = open(names_root + ".txt", 'w')
        names_file.write("\t".join(("ModelSEED ID","External ID","Source")) + "\n")
        for rxn in sorted(names_dict.keys()):
            for name in sorted(names_dict[rxn]):
                names_file.write("\t".join((rxn,name,'name')) + "\n")
        names_file.close()

    def saveAliases(self, alias_dict):
        alias_root = os.path.splitext(self.AliasFile)[0]

        # Print to TXT
        alias_file = open(alias_root + ".txt", 'w')
        alias_file.write("\t".join(("ModelSEED ID","External ID","Source")) + "\n")
        for rxn in sorted(alias_dict.keys()):
            for source in sorted (alias_dict[rxn].keys()):
                for alias in sorted(alias_dict[rxn][source]):
                    alias_file.write("\t".join((rxn,alias,source)) + "\n")
        alias_file.close()

    def saveReactions(self, reactions_dict):
        rxns_root = self.BiochemRoot + 'reaction_'

        # Initiate count
        rxns_count_thousands=0
        rxns_count_string=f"{rxns_count_thousands:02}"

        # Initiate TSV file handle
        rxns_tsv_file_handle = open(rxns_root + rxns_count_string+".tsv", 'w')
        rxns_tsv_file_handle.write("\t".join(self.Headers) + "\n")

        # Initiate JSON file handle
        rxns_json_file_handle = open(rxns_root + rxns_count_string+".json", 'w')
        rxns_json_list = list()

        # Initiate counting
        prev_rounded_count = 0
        
        # Iterate through reactions
        for rxn_id in sorted(reactions_dict.keys()):

            # Reset for every 1000
            # NB: we want a direct link between the filename and the reaction id
            # and for legacy reasons, there may not be a reaction id that falls on
            # a multiple of 1000, so, here we find the reaction id that "crosses"
            # a multiple of 1000 in order to keep count
            #
            # for the same legacy reasons, this means there won't be
            # 1000 reactions in every file
            rxn_count = int(rxn_id[3:])
            cur_rounded_count = rxn_count - rxn_count % 1000
            if(cur_rounded_count > prev_rounded_count):

                prev_rounded_count = cur_rounded_count

                # Write JSON list
                rxns_json_file_handle.write(json.dumps(rxns_json_list, indent=4, sort_keys=True))

                # Reset JSON list
                rxns_json_list = list()

                # Increment count
                rxns_count_thousands+=1
                rxns_count_string=f"{rxns_count_thousands:02}"

                # Reset tsv file handle
                rxns_tsv_file_handle.close()
                rxns_tsv_file_handle = open(rxns_root + rxns_count_string+".tsv", 'w')
                rxns_tsv_file_handle.write("\t".join(self.Headers) + "\n")

                # Reset json file handle
                rxns_json_file_handle.close()
                rxns_json_file_handle = open(rxns_root + rxns_count_string+".json", 'w')
                pass

            # Write TSV
            values_list=list()
            for header in self.Headers:
                value=reactions_dict[rxn_id][header]
                if(header=="stoichiometry"):
                    value = self.buildStoich(reactions_dict[rxn_id][header])
                elif(isinstance(value,list)):
                    value = "|".join(value)
                elif(isinstance(value,dict)):
                    entries = list()
                    for entry in value:
                        entries.append(entry+':'+value[entry])
                    value = "|".join(entries)
                values_list.append(str(value))
            rxns_tsv_file_handle.write("\t".join(values_list)+"\n")

            # Collect list for JSON
            rxn_obj = reactions_dict[rxn_id]
            for key in rxn_obj:
                if(isinstance(rxn_obj[key],dict)):
                    for entry in rxn_obj[key]:
                        if(rxn_obj[key][entry]=="null"):
                            rxn_obj[key][entry]=None
                if(rxn_obj[key]=="null"):
                    rxn_obj[key]=None

            rxns_json_list.append(rxn_obj)

        # Close TSV file handle
        rxns_tsv_file_handle.close()

        # Write last JSON list and close file handle
        rxns_json_file_handle.write(json.dumps(rxns_json_list, indent=4, sort_keys=True))
        rxns_json_file_handle.close()

    def loadMSAliases(self,sources_array=[]):
        if(len(sources_array)==0):
            sources_array.append("All")

        aliases_dict = dict()
        reader = DictReader(open(self.AliasFile), dialect = 'excel-tab')
        for line in reader:
            if("rxn" not in line['ModelSEED ID']):
                continue

            if("All" not in sources_array and line['Source'] not in sources_array):
                continue

            if(line['ModelSEED ID'] not in aliases_dict):
                   aliases_dict[line['ModelSEED ID']]=dict()

            for source in line['Source'].split('|'):
                if(source not in aliases_dict[line['ModelSEED ID']]):
                    aliases_dict[line['ModelSEED ID']][source]=list()

                aliases_dict[line['ModelSEED ID']][source].append(line['External ID'])

        return aliases_dict

    def loadNames(self):
        names_dict = dict()
        reader = DictReader(open(self.NameFile), dialect = 'excel-tab')
        for line in reader:
            if("rxn" not in line['ModelSEED ID']):
                continue

            if(line['ModelSEED ID'] not in names_dict):
                   names_dict[line['ModelSEED ID']]=list()

            names_dict[line['ModelSEED ID']].append(line['External ID'])

        return names_dict

    def loadPathways(self):
        pathways_dict = dict()
        reader = DictReader(open(self.PwyFile), dialect = 'excel-tab')
        for line in reader:
            if("rxn" not in line['ModelSEED ID']):
                continue

            if(line['ModelSEED ID'] not in pathways_dict):
                   pathways_dict[line['ModelSEED ID']]=dict()

            if(line['Source'] not in pathways_dict[line['ModelSEED ID']]):
                pathways_dict[line['ModelSEED ID']][line['Source']]=list()

            pathways_dict[line['ModelSEED ID']][line['Source']].append(line['External ID'])

        return pathways_dict

    def loadECs(self):
        ecs_dict = dict()
        reader = DictReader(open(self.ECFile), dialect = 'excel-tab')
        for line in reader:
            if("rxn" not in line['ModelSEED ID']):
                continue

            if(line['ModelSEED ID'] not in ecs_dict):
                   ecs_dict[line['ModelSEED ID']]=list()

            ecs_dict[line['ModelSEED ID']].append(line['External ID'])

        return ecs_dict
