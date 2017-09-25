import re, os, json

class Compounds:
    def __init__(self):
        self.BiochemRoot='../../Biochemistry/'
        self.CpdsFile=self.BiochemRoot+'compounds.tsv'

        CpdsFileHandle = open(self.CpdsFile,'r')
        HeaderLine = CpdsFileHandle.readline().strip()
        self.Headers = HeaderLine.split("\t")

    def loadCompounds(self):
        CpdsFile = open(self.CpdsFile,'r')
        CpdsFile.readline()

        CpdsDict = dict()
        for line in CpdsFile.readlines():
            line=line.strip()
            items = line.split("\t")
            CpdsDict[items[0]] = dict()
            for i in range(len(self.Headers)):
                item="null"
                if(i < len(items)):
                    item=items[i]

                #capture ints
                for header in ["is_core","is_obsolete","is_cofactor"]:
                    if(self.Headers[i] == header):
                        item=int(item)

                #capture floats
                #for header in ["mass","charge","deltag","deltagerr"]:
                #    if(self.Headers[i] == header):
                #        print header,item
                #        item=float(item)

                CpdsDict[items[0]][self.Headers[i]]=item

        return CpdsDict

    def parseFormula(self,Formula):
        if(Formula is None or Formula == "" or "noFormula" in Formula or "null" in Formula):
            return {} 
    
        Atoms = re.findall("\D[a-z]?\d*",Formula)

        Atoms_Dict = dict()
        for atom in Atoms:
            match = re.match("(\D[a-z]?)(\d*)",atom)
            Atoms_Dict[match.group(1)]=match.group(2)

            #Default empty string to 1
            if(Atoms_Dict[match.group(1)] == ""):
                Atoms_Dict[match.group(1)]=1
            else:
                Atoms_Dict[match.group(1)]=int(Atoms_Dict[match.group(1)])

        return Atoms_Dict

    def mergeFormula(self,Formula):
        Formula=Formula.strip()
        Notes = ""
        if(Formula is None or Formula == "" or "null" in Formula or len(re.findall("no[Ff]ormula",Formula))>0):
            return ("null",Notes)

        if(len(re.findall("(\)[nx])",Formula))>0):
            Notes="PO"
 
        Global_Atoms_Dict=dict()
        for subformula in re.findall("\(?([\w\s\.]+)\)?([nx*]?)?(\d?)",Formula):
            #The regex works, but returns empty hits for either beginning or end of string
            #The regex is trying to find formulas outside and within parentheses eg: Mg(Al,Fe)Si4O10(OH).4H2O
            subformula_string = subformula[0].strip()
            if(subformula_string != ''):
                bracketed_multiplier = 1
                #Redundant but worth being explicit: generic polymeric formulas assumed to be 1 unit
                if(len(re.findall("[nx*]",subformula[1]))==0 and subformula[2] != ""):
                    bracketed_multiplier = int(subformula[2])

                #Avoid empty strings
                for fragment in (x for x in subformula_string.split(".") if x):
                    fragment=fragment.strip()
                    fragment_multiplier=1
                    #Fragments can have a multiplier at the beginning of the string, such as 4H2O
                    if(len(re.findall("^(\d+)(.*)$",fragment))):
                        (fragment_multiplier,fragment)=re.findall("^(\d+)(.*)$",fragment)[0]
                        fragment_multiplier=int(fragment_multiplier)

                    Fragment_Atoms_Dict = self.parseFormula(fragment)
                    for atom in Fragment_Atoms_Dict:
                        if atom not in Global_Atoms_Dict.keys():
                            Global_Atoms_Dict[atom]=0
                        Global_Atoms_Dict[atom]+=Fragment_Atoms_Dict[atom]*bracketed_multiplier*fragment_multiplier

        return (self.buildFormula(Global_Atoms_Dict),Notes)

    def buildFormula(self,Atoms_Dict):
        Formula=""
        for atom in self.hill_sorted(Atoms_Dict.keys()):
            if(Atoms_Dict[atom]==1):
                Atoms_Dict[atom]=""
            Formula+=atom+str(Atoms_Dict[atom])
        return Formula

    def hill_sorted(self,atoms):
        if("C" in atoms):
            atoms.remove("C")
            yield "C"
        if("H" in atoms):
            atoms.remove("H")
            yield "H"
        for atom in sorted(atoms):
            yield atom

    def saveCompounds(self,Compounds_Dict):
        CpdsRoot = os.path.splitext(self.CpdsFile)[0]
        
        #Print to TSV
        CpdsFile = open(CpdsRoot+".tsv",'w')
        CpdsFile.write("\t".join(self.Headers)+"\n")
        for cpd in sorted(Compounds_Dict.keys()):
            CpdsFile.write("\t".join(str(Compounds_Dict[cpd][header]) for header in self.Headers)+"\n")
        CpdsFile.close()

        #Print to JSON
        CpdsFile = open(CpdsRoot+".json",'w')
        CpdsFile.write(json.dumps(Compounds_Dict, indent=4, sort_keys=True))
        CpdsFile.close()
