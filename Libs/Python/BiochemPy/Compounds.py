import re
import os
import json
import glob
from csv import DictReader

class Compounds:
    def __init__(self, biochem_root='../../../Biochemistry/',
                 cpds_file='compound_00.tsv'):

        self.BiochemRoot = os.path.dirname(__file__)+'/'+biochem_root
        self.CpdsFile = self.BiochemRoot + cpds_file
        self.AliasFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Compound_Aliases.txt"
        self.NameFile = self.BiochemRoot + "Aliases/Unique_ModelSEED_Compound_Names.txt"
        self.StructRoot = self.BiochemRoot + "Structures/"

        reader = DictReader(open(self.CpdsFile), dialect='excel-tab')
        self.Headers = reader.fieldnames

    def loadCompounds(self):

        search_path = os.path.join(self.BiochemRoot,"compound_*.json")
        cpds_dict = dict()
        for compounds_file in sorted(glob.glob(search_path)):
            with open(compounds_file) as json_file_handle:
                cpds_list = json.load(json_file_handle)
                for cpd_obj in cpds_list:
                    for key in cpd_obj:
                        if(isinstance(cpd_obj[key],list)):
                            for i in range(len(cpd_obj[key])):
                                if(cpd_obj[key][i] is None):
                                    cpd_obj[key][i]="null"
                            
                        if(isinstance(cpd_obj[key],dict)):
                            for entry in cpd_obj[key]:
                                if(cpd_obj[key][entry] is None):
                                    cpd_obj[key][entry]="null"
                            
                        if(cpd_obj[key] is None):
                            cpd_obj[key]="null"
                        
                    cpds_dict[cpd_obj['id']]=cpd_obj
                    
        return cpds_dict

    def loadCompounds_tsv(self):
        print("WARNING: This function is currently redundant and will only load one file!")
        reader = DictReader(open(self.CpdsFile), dialect='excel-tab')
        type_mapping = {"is_core": int, "is_obsolete": int, "is_cofactor": int, "charge": int,
                        "mass": float, "deltag": float, "deltagerr": float}
        lists = ["aliases","notes"]
        dicts = []

        cpds_dict = dict()
        for line in reader:
            for list_type in lists:
                if(line[list_type] != "null"):
                    if(line[list_type] == ""):
                        line[list_type]="null"
                    else:
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
            fields_array= ['ID','Source','Alias','Formula','Charge','Structure']

            if(unique==False):
                struct_file = "All_ModelSEED_Structures.txt"
                fields_array = ['ID','Source','Type','Alias','DB','Formula','Charge','Structure']

            struct_file = self.StructRoot+struct_file
            reader = DictReader(open(struct_file), dialect = "excel-tab", fieldnames = fields_array)
            for line in reader:
                if("cpd" not in line['ID']):
                    continue

                if(line['ID'] not in structures_dict):
                    structures_dict[line['ID']]={}

                if(line['Source'] in sources_array):
                    if(line['Source'] not in structures_dict[line['ID']]):
                        structures_dict[line['ID']][line['Source']]=dict()
                    structures_dict[line['ID']][line['Source']][line['Structure']]={'formula':line['Formula'],
                                                                                    'charge':line['Charge'],
                                                                                    'alias':line['Alias'].split(';')}
                    if('Type' in line):
                        structures_dict[line['ID']][line['Source']][line['Structure']]['type']=line['Type']
            return structures_dict

        # Per-source layout (post-A1 migration):
        #   <db>/inchi.tsv, smiles.tsv, inchikey.tsv  → "Original" stage
        #   <db>/protonations/*.tsv                   → "Charged" stage
        # Stage names "Original" / "Charged" preserved so callers
        # (List_ModelSEED_Structures.py, etc.) keep working unchanged.
        original_files   = {'InChI': 'inchi.tsv', 'SMILE': 'smiles.tsv', 'InChIKey': 'inchikey.tsv'}
        original_columns = {'InChI': 'inchi',     'SMILE': 'smiles',     'InChIKey': 'inchikey'}

        for struct_type in sources_array:
            structures_dict[struct_type] = dict()
            for db in db_array:
                # ---- "Original" stage: source-as-downloaded ----
                f = self.StructRoot + db + "/" + original_files.get(struct_type, '')
                if os.path.isfile(f):
                    with open(f) as fh:
                        reader = DictReader(fh, dialect='excel-tab')
                        col = original_columns[struct_type]
                        for line in reader:
                            ext_id = line.get('external_id') or line.get('ID')
                            struct = line.get(col)
                            if not ext_id or not struct:
                                continue
                            structures_dict[struct_type].setdefault(ext_id, {}).setdefault('Original', {})[struct] = 1

                # ---- "Charged" stage: protonations/<tool>_<ver>_ph<n>.tsv ----
                proto_dir = self.StructRoot + db + "/protonations"
                if os.path.isdir(proto_dir):
                    for proto_file in sorted(glob.glob(proto_dir + "/*.tsv")):
                        with open(proto_file) as fh:
                            reader = DictReader(fh, dialect='excel-tab')
                            for line in reader:
                                if line.get('type') != struct_type:
                                    continue
                                ext_id = line.get('external_id')
                                struct = line.get('structure')
                                if not ext_id or not struct:
                                    continue
                                structures_dict[struct_type].setdefault(ext_id, {}).setdefault('Charged', {})[struct] = 1

        return structures_dict

    def loadPerSourceFormulasCharges(self, struct_types=None, db_array=None):
        """Returns formulas[db][struct_type][stage][ext_id] = {'formula', 'charge'}.

        Reads from the new layout: inchi.tsv / smiles.tsv carry per-source
        Original-stage formula+charge; protonations/*.tsv carries Charged-stage
        formula+charge filtered by the 'type' column. Mirrors the dict shape
        that List_ModelSEED_Structures.py used to build by reading
        *_Formulas_Charges.txt files directly.
        """
        if struct_types is None:
            struct_types = ['InChI', 'SMILE']
        if db_array is None:
            db_array = ['KEGG', 'MetaCyc', 'ChEBI', 'Rhea']

        original_files   = {'InChI': 'inchi.tsv', 'SMILE': 'smiles.tsv'}
        out = {}
        for db in db_array:
            out[db] = {}
            for struct_type in struct_types:
                out[db][struct_type] = {'Charged': {}, 'Original': {}}

                # Original
                if struct_type in original_files:
                    f = self.StructRoot + db + '/' + original_files[struct_type]
                    if os.path.isfile(f):
                        with open(f) as fh:
                            reader = DictReader(fh, dialect='excel-tab')
                            for line in reader:
                                ext_id = line.get('external_id')
                                if not ext_id:
                                    continue
                                if line.get('formula') or line.get('charge'):
                                    out[db][struct_type]['Original'][ext_id] = {
                                        'formula': line.get('formula', ''),
                                        'charge':  line.get('charge', ''),
                                    }

                # Charged
                proto_dir = self.StructRoot + db + '/protonations'
                if os.path.isdir(proto_dir):
                    for proto_file in sorted(glob.glob(proto_dir + '/*.tsv')):
                        with open(proto_file) as fh:
                            reader = DictReader(fh, dialect='excel-tab')
                            for line in reader:
                                if line.get('type') != struct_type:
                                    continue
                                ext_id = line.get('external_id')
                                if not ext_id:
                                    continue
                                if line.get('formula') or line.get('charge'):
                                    out[db][struct_type]['Charged'][ext_id] = {
                                        'formula': line.get('formula', ''),
                                        'charge':  line.get('charge', ''),
                                    }
        return out

    def loadPerSourceCharges(self, db_array=None, ph=7):
        """Returns charges[(db, ext_id)] = int(charge_at_ph).

        Reads from <db>/protonations/marvin_*_ph<n>.tsv (the post-A1 layout).
        Iterates every matching TSV in the protonations/ directory; if a source
        ships multiple snapshots, the last-sorted one wins for each ext_id.

        Each compound typically has multiple rows in the protonations file
        (one per structure type: InChI, SMILE, InChIKey) but they all carry
        the same charge — last-write-wins collapses them safely.

        Applies the same ID normalization as loadPerSourcePkas so the
        returned ext_id matches the alias-file convention.
        """
        if db_array is None:
            db_array = ['KEGG', 'MetaCyc', 'ChEBI', 'Rhea']

        out = {}
        for db in db_array:
            proton_dir = self.StructRoot + db + '/protonations'
            if not os.path.isdir(proton_dir):
                continue
            pattern = proton_dir + '/*_ph' + str(ph) + '.tsv'
            for proton_file in sorted(glob.glob(pattern)):
                with open(proton_file) as fh:
                    reader = DictReader(fh, dialect='excel-tab')
                    for line in reader:
                        ext_id = line.get('external_id')
                        charge = line.get('charge')
                        if not ext_id or charge in (None, ''):
                            continue
                        if db == 'ChEBI' and ext_id.startswith('CHEBI_'):
                            ext_id = ext_id[len('CHEBI_'):]
                        elif db == 'Rhea' and ext_id.startswith('POLYMER_'):
                            ext_id = 'POLYMER:' + ext_id[len('POLYMER_'):]
                        try:
                            out[(db, ext_id)] = int(charge)
                        except ValueError:
                            continue
        return out

    def loadPerSourcePkas(self, db_array=None):
        """Returns pkas[(db, ext_id)] = {'pKa': value, 'pKb': value}.

        Reads from <db>/pkas/*.tsv (the post-A1 layout). Iterates every TSV
        in the pkas/ directory; if a source ships multiple snapshots, the
        last-sorted one wins for each (ext_id, kind) pair.

        Normalizes source-specific id quirks so that the returned ext_id
        matches the alias-file convention:
          - ChEBI pKa rows are keyed 'CHEBI_15377' but Aliases lists '15377'.
          - Rhea pKa rows are keyed 'POLYMER_10033' but Aliases lists
            'POLYMER:10033'.
        Without this normalization the lookup never matches and the
        pKa data is silently unused. See sources.yaml for the
        consumed_by_production flag history.
        """
        if db_array is None:
            db_array = ['KEGG', 'MetaCyc', 'ChEBI', 'Rhea']

        out = {}
        for db in db_array:
            pka_dir = self.StructRoot + db + '/pkas'
            if not os.path.isdir(pka_dir):
                continue
            for pka_file in sorted(glob.glob(pka_dir + '/*.tsv')):
                with open(pka_file) as fh:
                    reader = DictReader(fh, dialect='excel-tab')
                    for line in reader:
                        ext_id = line.get('external_id')
                        kind   = line.get('kind')
                        value  = line.get('value')
                        if not ext_id or not kind:
                            continue
                        if db == 'ChEBI' and ext_id.startswith('CHEBI_'):
                            ext_id = ext_id[len('CHEBI_'):]
                        elif db == 'Rhea' and ext_id.startswith('POLYMER_'):
                            ext_id = 'POLYMER:' + ext_id[len('POLYMER_'):]
                        out.setdefault((db, ext_id), {})[kind] = value
        return out

    @staticmethod
    def searchname(name):
        searchnames_list = [name]
        searchname = name.lower()
        searchnames_list.append(searchname)
        
        #try to keep/maintain charges
        ending = ""
        if(searchname.endswith("-")):
            ending="-"

        if(searchname.endswith("+")):
            ending="+"

        searchname = ''.join(char for char in searchname if char.isalnum())
        searchnames_list.append(searchname+ending)
        
        #attempting to match fatty acids
        if(re.search('icacid$',searchname)):
            searchname = re.sub('icacid','ate',searchname)
            searchnames_list.append(searchname+ending)
        elif(re.search('ate$',searchname)):
            searchname = re.sub('ate','icacid',searchname)
            searchnames_list.append(searchname+ending)
            
        #remove redundant articles
        if(re.search(r'^an?\s',searchname)):
            searchname = re.sub('^an?','',searchname)
            searchnames_list.append(searchname+ending)

        return searchnames_list

    @staticmethod
    def parseFormula(formula):
        if (formula.strip() in {None, "", "noFormula", "null"}):
            return {}

        atoms = re.findall(r"\D[a-z]?\d*", formula)

        atoms_dict = dict()
        for atom in atoms:
            match = re.match(r"(\D[a-z]?)(\d*)", atom)
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

        if (len(re.findall(r"(\)[nx])", formula)) > 0):
            Notes = "PO"

        global_atoms_dict = dict()
        for subformula in re.findall(r"\(?([\w\s\.]+)\)?([nx*]?)?(\d?)",
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
                    if (len(re.findall(r"^(\d+)(.*)$", fragment))):
                        (fragment_multiplier, fragment) = \
                        re.findall(r"^(\d+)(.*)$", fragment)[0]
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

        cpds_root = self.BiochemRoot + 'compound_'

        # Initiate count
        cpds_count_thousands=0
        cpds_count_string=f"{cpds_count_thousands:02}"

        # Initiate TSV file handle
        cpds_tsv_file_handle = open(cpds_root + cpds_count_string+".tsv", 'w')
        cpds_tsv_file_handle.write("\t".join(self.Headers) + "\n")

        # Initiate JSON file handle
        cpds_json_file_handle = open(cpds_root + cpds_count_string+".json", 'w')
        cpds_json_list = list()

        # Initiate counting
        prev_rounded_count = 0
        
        # Iterate through compounds
        for cpd_id in sorted(compounds_dict.keys()):

            # Reset for every 1000
            # NB: we want a direct link between the filename and the compound id
            # and for legacy reasons, there may not be a compound id that falls on
            # a multiple of 1000, so, here we find the compound id that "crosses"
            # a multiple of 1000 in order to keep count
            #
            # for the same legacy reasons, this means there won't be
            # 1000 compounds in every file
            cpd_count = int(cpd_id[3:])
            cur_rounded_count = cpd_count - cpd_count % 1000
            if(cur_rounded_count > prev_rounded_count):

                prev_rounded_count = cur_rounded_count

                # Write JSON list
                cpds_json_file_handle.write(json.dumps(cpds_json_list, indent=4, sort_keys=True))

                # Reset JSON list
                cpds_json_list = list()

                # Increment count
                cpds_count_thousands+=1
                cpds_count_string=f"{cpds_count_thousands:02}"

                # Reset tsv file handle
                cpds_tsv_file_handle.close()
                cpds_tsv_file_handle = open(cpds_root + cpds_count_string+".tsv", 'w')
                cpds_tsv_file_handle.write("\t".join(self.Headers) + "\n")

                # Reset json file handle
                cpds_json_file_handle.close()
                cpds_json_file_handle = open(cpds_root + cpds_count_string+".json", 'w')
                pass

            # Write TSV
            values_list=list()
            for header in self.Headers:
                value=compounds_dict[cpd_id][header]
                if(value is None):
                    value="null"
                if(isinstance(value,list)):
                    value = "|".join(value)
                if(isinstance(value,dict)):
                    entries = list()
                    for entry in value:
                        entries.append(entry+':'+value[entry])
                    value = "|".join(entries)
                values_list.append(str(value))
            cpds_tsv_file_handle.write("\t".join(values_list)+"\n")

            # Collect list for JSON
            cpd_obj = compounds_dict[cpd_id]
            for key in cpd_obj:
                if(isinstance(cpd_obj[key],dict)):
                    for entry in cpd_obj[key]:
                        if(cpd_obj[key][entry]=="null"):
                            cpd_obj[key][entry]=None
                if(cpd_obj[key]=="null"):
                    cpd_obj[key]=None

            cpds_json_list.append(cpd_obj)

        # Close TSV file handle
        cpds_tsv_file_handle.close()

        # Write last JSON list and close file handle
        cpds_json_file_handle.write(json.dumps(cpds_json_list, indent=4, sort_keys=True))
        cpds_json_file_handle.close()
