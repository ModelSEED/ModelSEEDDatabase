import os
import json
from csv import DictReader


class Reactions:
    def __init__(self, biochem_root='../../Biochemistry/',
                 rxns_file='reactions.tsv'):
        self.BiochemRoot = biochem_root
        self.RxnsFile = biochem_root + rxns_file
        self.AliasFile = biochem_root + "Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
        self.NameFile = biochem_root + "Aliases/Unique_ModelSEED_Reaction_Names.txt"
        self.ECFile = biochem_root + "Aliases/Unique_ModelSEED_Reaction_ECs.txt"

        reader = DictReader(open(self.RxnsFile), dialect='excel-tab')
        self.Headers = reader.fieldnames

        from BiochemPy import Compounds
        self.CompoundsHelper = Compounds()
        self.Compounds_Dict = self.CompoundsHelper.loadCompounds()

    def loadReactions(self):
        reader = DictReader(open(self.RxnsFile), dialect='excel-tab')
        rxns_dict = dict()
        for line in reader:
            for header in ["is_transport", "is_obsolete"]:
                line[header] = int(line[header])
            rxns_dict[line['id']] = line

        return rxns_dict

    def parseStoich(self, stoichiometry):
        rxn_cpds_array = list()
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

    @staticmethod
    def isTransport(rxn_cpds_array):
        compartments_dict=dict()
        for rgt in rxn_cpds_array:
            compartments_dict[rgt['compartment']]=1
        if(len(compartments_dict.keys())>1):
            return 1
        else:
            return 0

    def generateCodes(self, rxns_dict):
        codes_dict=dict()
        for rxn in rxns_dict:
            if(rxns_dict[rxn]['status']=="EMPTY"):
                continue
            code = self.generateCode(rxns_dict[rxn]['stoichiometry'])
            if(code not in codes_dict):
                codes_dict[code]=dict()
            codes_dict[code][rxn]=1
        return codes_dict

    def generateCode(self,stoichiometry):
        rxn_cpds_array = self.parseStoich(stoichiometry)

        #It matters if its a transport reaction, and we include protons when matching transpor
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

            rgt["coefficient"] = str(rgt["coefficient"])
            rgt["compartment"] = str(rgt["compartment"])
            rgt["index"] = str(rgt["index"])

            rgt_string = ":".join(
                [rgt["coefficient"], rgt["compound"], rgt["compartment"],
                 rgt["index"], rgt["name"]])
            stoichiometry_array.append(rgt_string)
        stoichiometry_string = ";".join(stoichiometry_array)
        return stoichiometry_string

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
                return "ERROR: Duplicate reagents"

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
            rgt["coefficient"] = cpds_coeff_dict[rgt["compound"]]
            cpds_dict[rgt["compound"]] = rgt

        ########################################
        # Check for duplicate elements, across
        # all compounds, these are balanced 
        # directly.
        #######################################
        rxn_net_charge = 0.0
        rxn_net_mass = dict()
        for cpd in cpds_dict.keys():
            cpd_atoms = self.CompoundsHelper.parseFormula(
                cpds_dict[cpd]["formula"])

            if (len(cpd_atoms.keys()) == 0):
                return "CPDFORMERROR"

            cpd_coeff_charge = float(cpds_dict[cpd]["charge"]) * float(
                cpds_dict[cpd]["coefficient"])
            rxn_net_charge += cpd_coeff_charge

            for atom in cpd_atoms.keys():
                atom_coeff_mass = float(cpd_atoms[atom]) * float(
                    cpds_dict[cpd]["coefficient"])

                if (atom not in rxn_net_mass.keys()):
                    rxn_net_mass[atom] = 0.0

                rxn_net_mass[atom] += atom_coeff_mass

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

        return

    def rebuildReaction(self, reaction_dict, stoichiometry):
        # Assign stoich
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

        return

    def saveReactions(self, reactions_dict):
        rxns_root = os.path.splitext(self.RxnsFile)[0]

        # Print to TSV
        rxns_file = open(rxns_root + ".tsv", 'w')
        rxns_file.write("\t".join(self.Headers) + "\n")
        for rxn in sorted(reactions_dict.keys()):
            rxns_file.write("\t".join(
                str(reactions_dict[rxn][header]) for header in
                self.Headers) + "\n")
        rxns_file.close()

        # Print to JSON
        rxns_file = open(rxns_root + ".json", 'w')
        rxns_file.write(json.dumps(reactions_dict, indent=4, sort_keys=True))
        rxns_file.close()

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
                   names_dict[line['ModelSEED ID']]=dict()

            #redundant as only one source but keep this just in case
            for source in line['Source'].split('|'):
                if(source not in names_dict[line['ModelSEED ID']]):
                    names_dict[line['ModelSEED ID']][source]=list()

                names_dict[line['ModelSEED ID']][source].append(line['External ID'])

        return names_dict

    def loadECs(self):
        ecs_dict = dict()
        reader = DictReader(open(self.ECFile), dialect = 'excel-tab')
        for line in reader:
            if("rxn" not in line['ModelSEED ID']):
                continue

            if(line['ModelSEED ID'] not in ecs_dict):
                   ecs_dict[line['ModelSEED ID']]=dict()

            #redundant as only one source but keep this just in case
            for source in line['Source'].split('|'):
                if(source not in ecs_dict[line['ModelSEED ID']]):
                    ecs_dict[line['ModelSEED ID']][source]=list()

                ecs_dict[line['ModelSEED ID']][source].append(line['External ID'])

        return ecs_dict
