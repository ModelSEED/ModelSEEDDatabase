import re
import os
import json
from csv import DictReader

class Compounds:
    def __init__(self, biochem_root='../../../Biochemistry/',
                 cpds_file='compounds.tsv'):

        self.BiochemRoot = os.path.dirname(__file__)+'/'+biochem_root
        self.CpdsFile = self.BiochemRoot + cpds_file
        self.AliasFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Compound_Aliases.txt"
        self.NameFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Compound_Names.txt"
        self.StructRoot = self.BiochemRoot + "Structures/"

        reader = DictReader(open(self.CpdsFile), dialect='excel-tab')
        self.Headers = reader.fieldnames

    def loadCompounds(self):
        reader = DictReader(open(self.CpdsFile), dialect='excel-tab')
        cpds_dict = {}
        for line in reader:
            for header in ["is_core", "is_obsolete", "is_cofactor"]:
                line[header] = int(line[header])
            cpds_dict[line['id']] = line

        return cpds_dict

    def loadMSAliases(self,sources_array=[]):
        if(len(sources_array)==0):
            sources_array.append("All")

        aliases_dict = dict()
        reader = DictReader(open(self.AliasFile), dialect = 'excel-tab')
        for line in reader:
            if("cpd" not in line['ModelSEED ID']):
                continue

            for source in line['Source'].split('|'):

                if("All" not in sources_array and source not in sources_array):
                    continue

                if(line['ModelSEED ID'] not in aliases_dict):
                   aliases_dict[line['ModelSEED ID']]=dict()

                if(source not in aliases_dict[line['ModelSEED ID']]):
                    aliases_dict[line['ModelSEED ID']][source]=list()

                aliases_dict[line['ModelSEED ID']][source].append(line['External ID'])

        return aliases_dict

    def loadSourceAliases(self):
        aliases_dict = dict()
        reader = DictReader(open(self.AliasFile), dialect = 'excel-tab')
        for line in reader:
            if("cpd" not in line['ModelSEED ID']):
                continue

            for source in line['Source'].split('|'):
                if(source not in aliases_dict):
                    aliases_dict[source]=dict()

                if(line['External ID'] not in aliases_dict[source]):
                    aliases_dict[source][line['External ID']]=list()

                aliases_dict[source][line['External ID']].append(line['ModelSEED ID'])

        return aliases_dict

    def loadNames(self):
        names_dict = dict()
        reader = DictReader(open(self.NameFile), dialect = 'excel-tab')
        for line in reader:
            if("cpd" not in line['ModelSEED ID']):
                continue

            if(line['ModelSEED ID'] not in names_dict):
                   names_dict[line['ModelSEED ID']]=list()

            names_dict[line['ModelSEED ID']].append(line['External ID'])

        return names_dict

    def loadStructures(self,sources_array=[],db_array=[],unique=True):
        if(len(sources_array)==0):
            sources_array=["SMILE","InChIKey","InChI"]

        if(len(db_array)==0):
            db_array=["KEGG","MetaCyc"]

        structures_dict = dict()
        if(len(db_array)==1 and db_array[0]=="ModelSEED"):
            struct_file = "Unique_ModelSEED_Structures.txt"
            fields_array= ['ID','Source','Aliases','Formula','Charge','Structure']

            if(unique==False):
                struct_file = "All_ModelSEED_Structures.txt"
                fields_array = ['ID','Source','Type','Alias','DB','Formula','Charge','Structure']

            struct_file = self.StructRoot+struct_file
            reader = DictReader(open(struct_file), dialect = "excel-tab", fieldnames = fields_array)
            for line in reader:
                if(line['ID'] not in structures_dict):
                    structures_dict[line['ID']]={}

                if(line['Source'] in sources_array):
                    if(line['Source'] not in structures_dict[line['ID']]):
                        structures_dict[line['ID']][line['Source']]=dict()
                    structures_dict[line['ID']][line['Source']][line['Structure']]={'formula':line['Formula'],
                                                                                    'charge':line['Charge']}

            return structures_dict

        for struct_type in sources_array:
            structures_dict[struct_type]=dict()
            for db in db_array:
                for struct_stage in ["Charged","Original"]:
                    struct_file = db+"/"+struct_type+"_"+struct_stage+"Strings.txt"
                    struct_file = self.StructRoot+struct_file

                    if(os.path.isfile(struct_file)==False):
                        continue

                    reader = DictReader(open(struct_file), dialect = "excel-tab", fieldnames = ['ID','Structure','Name'])
                    for line in reader:
                        if(line['ID'] not in structures_dict[struct_type]):
                            structures_dict[struct_type][line['ID']]=dict()

                        if(struct_stage not in structures_dict[struct_type][line['ID']]):
                            structures_dict[struct_type][line['ID']][struct_stage]=dict()

                        structures_dict[struct_type][line['ID']][struct_stage][line['Structure']]=1

        return structures_dict

    @staticmethod
    def searchname(name):
        searchname = name.lower()

        #try to keep/maintain charges
        ending = ""
        if(searchname.endswith("-")):
            ending="-"

        if(searchname.endswith("+")):
            ending="+"

        searchname = ''.join(char for char in searchname if char.isalnum())

        #attempting to match fatty acids
        searchname = re.sub('icacid','ate',searchname)

        #remove redundant articles
        if(re.search('^an?\s',name)):
            searchname = re.sub('^an?','',searchname)

        searchname+=ending

        return searchname

    @staticmethod
    def parseFormula(formula):
        if (formula.strip() in {None, "", "noFormula", "null"}):
            return {}

        atoms = re.findall("\D[a-z]?\d*", formula)

        atoms_dict = dict()
        for atom in atoms:
            match = re.match("(\D[a-z]?)(\d*)", atom)
            atoms_dict[match.group(1)] = match.group(2)

            # Default empty string to 1
            if (atoms_dict[match.group(1)] == ""):
                atoms_dict[match.group(1)] = 1
            else:
                atoms_dict[match.group(1)] = int(atoms_dict[match.group(1)])

        return atoms_dict

    @staticmethod
    def mergeFormula(formula):
        formula = formula.strip()
        Notes = ""
        if (formula is None or formula == "" or "null" in formula or len(
                re.findall("no[Ff]ormula", formula)) > 0):
            return ("null", Notes)

        if (len(re.findall("(\)[nx])", formula)) > 0):
            Notes = "PO"

        global_atoms_dict = dict()
        for subformula in re.findall("\(?([\w\s\.]+)\)?([nx*]?)?(\d?)",
                                     formula):
            # The regex works, but returns empty hits for either beginning or end of string
            # The regex is trying to find formulas outside and within parentheses eg: Mg(Al,Fe)Si4O10(OH).4H2O
            subformula_string = subformula[0].strip()
            if (subformula_string != ''):
                bracketed_multiplier = 1
                # Redundant but worth being explicit: generic polymeric formulas assumed to be 1 unit
                if (len(re.findall("[nx*]", subformula[1])) == 0 and
                            subformula[2] != ""):
                    bracketed_multiplier = int(subformula[2])

                # Avoid empty strings
                for fragment in (x for x in subformula_string.split(".") if x):
                    fragment = fragment.strip()
                    fragment_multiplier = 1
                    # Fragments can have a multiplier at the beginning of the string, such as 4H2O
                    if (len(re.findall("^(\d+)(.*)$", fragment))):
                        (fragment_multiplier, fragment) = \
                        re.findall("^(\d+)(.*)$", fragment)[0]
                        fragment_multiplier = int(fragment_multiplier)

                    fragment_atoms_dict = Compounds.parseFormula(fragment)
                    for atom in fragment_atoms_dict:
                        if atom not in global_atoms_dict.keys():
                            global_atoms_dict[atom] = 0
                        global_atoms_dict[atom] += fragment_atoms_dict[atom] \
                            * bracketed_multiplier * fragment_multiplier

        return (Compounds.buildFormula(global_atoms_dict), Notes)

    @staticmethod
    def buildFormula(Atoms_Dict):
        formula = ""
        for atom in Compounds.hill_sorted(list(Atoms_Dict.keys())):
            if (Atoms_Dict[atom] == 1):
                Atoms_Dict[atom] = ""
            formula += atom + str(Atoms_Dict[atom])
        return formula

    @staticmethod
    def hill_sorted(atoms):
        if ("C" in atoms):
            atoms.remove("C")
            yield "C"
        if ("H" in atoms):
            atoms.remove("H")
            yield "H"
        for atom in sorted(atoms):
            yield atom

    def saveNames(self, names_dict):
        names_root = os.path.splitext(self.NameFile)[0]

        # Print to TXT
        names_file = open(names_root + ".txt", 'w')
        names_file.write("\t".join(("ModelSEED ID","External ID","Source")) + "\n")
        for cpd in sorted(names_dict.keys()):
            for name in sorted(names_dict[cpd]):
                names_file.write("\t".join((cpd,name,'name')) + "\n")
        names_file.close()

    def saveAliases(self, alias_dict):
        alias_root = os.path.splitext(self.AliasFile)[0]

        # Print to TXT
        alias_file = open(alias_root + ".txt", 'w')
        alias_file.write("\t".join(("ModelSEED ID","External ID","Source")) + "\n")
        for cpd in sorted(alias_dict.keys()):
            for source in sorted (alias_dict[cpd].keys()):
                for alias in sorted(alias_dict[cpd][source]):
                    alias_file.write("\t".join((cpd,alias,source)) + "\n")
        alias_file.close()

    def saveCompounds(self, compounds_dict):
        cpds_root = os.path.splitext(self.CpdsFile)[0]

        # Print to TSV
        cpds_file = open(cpds_root + ".tsv", 'w')
        cpds_file.write("\t".join(self.Headers) + "\n")
        for cpd in sorted(compounds_dict.keys()):
            cpds_file.write("\t".join(
                str(compounds_dict[cpd][header]) for header in
                self.Headers) + "\n")
        cpds_file.close()

        #Re-configure JSON
        new_compounds_dict = list()
        for cpd_id in sorted(compounds_dict):
            cpd_obj = compounds_dict[cpd_id]
            for key in cpd_obj:
                if(cpd_obj[key]=="null"):
                    cpd_obj[key]=None
            new_compounds_dict.append(cpd_obj)

        # Print to JSON
        cpds_file = open(cpds_root + ".json", 'w', newline='\n')
        cpds_file.write(json.dumps(new_compounds_dict, indent=4, sort_keys=True))
        cpds_file.close()
