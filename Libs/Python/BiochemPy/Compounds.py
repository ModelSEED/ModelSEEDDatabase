import re
import os
import json
from csv import DictReader


class Compounds:
    def __init__(self, biochem_root='../../Biochemistry/',
                 cpds_file='compounds.tsv'):
        self.BiochemRoot = biochem_root
        self.CpdsFile = biochem_root + cpds_file

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

    @staticmethod
    def parseFormula(formula):
        if (
                        formula is None or formula == "" or "noFormula" in formula or "null" in formula):
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
        for atom in Compounds.hill_sorted(Atoms_Dict.keys()):
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

        # Print to JSON
        cpds_file = open(cpds_root + ".json", 'w')
        cpds_file.write(json.dumps(compounds_dict, indent=4, sort_keys=True))
        cpds_file.close()
