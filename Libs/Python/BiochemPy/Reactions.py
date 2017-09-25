import os, json

class Reactions:
    def __init__(self):
        self.BiochemRoot='../../Biochemistry/'
        self.RxnsFile=self.BiochemRoot+'reactions.tsv'

        RxnsFileHandle = open(self.RxnsFile,'r')
        HeaderLine = RxnsFileHandle.readline().strip()
        self.Headers = HeaderLine.split("\t")

        from BiochemPy import Compounds
        self.CompoundsHelper=Compounds()
        self.Compounds_Dict = self.CompoundsHelper.loadCompounds()

    def loadReactions(self):
        RxnsFile = open(self.RxnsFile,'r')
        RxnsFile.readline()

        RxnsDict = dict()
        for line in RxnsFile.readlines():
            line=line.strip()
            items = line.split("\t")
            RxnsDict[items[0]] = dict()
            for i in range(len(self.Headers)):
                item="null"
                if(i < len(items)):
                    item=items[i]

                #capture ints
                for header in ["is_transport","is_obsolete"]:
                    if(self.Headers[i] == header):
                        item=int(item)

                #capture floats
                #for header in ["deltag","deltagerr"]:
                #    if(self.Headers[i] == header):
                #        item=float(item)

                RxnsDict[items[0]][self.Headers[i]]=item

        return RxnsDict

    def parseStoich(self,stoichiometry):
        Rxn_Cpds_Array=list()
        for rgt in stoichiometry.split(";"):
            (coeff,cpd,cpt,index,name) = rgt.split(":",4)
            rgt_id = cpd+"_"+cpt+index

            coeff=float(coeff)

            #Correct for redundant ".0" in floats
            if(str(coeff)[-2:] == ".0"):
                coeff=int(round(coeff))

            cpt=int(cpt)
            index=int(index)

            Rxn_Cpds_Array.append({"reagent":rgt_id,"coefficient":coeff,
                                   "compound":cpd,"compartment":cpt,"index":index,"name":name,
                                   "formula":self.Compounds_Dict[cpd]["formula"],
                                   "charge":self.Compounds_Dict[cpd]["charge"]})
        return Rxn_Cpds_Array

    def buildStoich(self,Rxn_Cpds_Array):
        Stoichiometry_Array = list()
        for rgt in sorted(Rxn_Cpds_Array, key = lambda x: (int(x["coefficient"] > 0),x["reagent"])):

            #Correct for redundant ".0" in floats
            if(str(rgt["coefficient"])[-2:] == ".0"):
                rgt["coefficient"]=int(round(rgt["coefficient"]))

            rgt["coefficient"]=str(rgt["coefficient"])
            rgt["compartment"]=str(rgt["compartment"])
            rgt["index"]=str(rgt["index"])

            Rgt_String = ":".join([rgt["coefficient"],rgt["compound"],rgt["compartment"],rgt["index"],rgt["name"]])
            Stoichiometry_Array.append(Rgt_String)
        Stoichiometry_String = ";".join(Stoichiometry_Array)
        return Stoichiometry_String

    def balanceReaction(self,Rgts_Array):
        if(len(Rgts_Array)==0):
            return "EMPTY"

        ########################################
        # Check that each reagent is either a 
        # different compound or in a different
        # compartment, and report.
        ########################################
        Rgts_Dict=dict()
        for rgt in Rgts_Array:
            if(rgt["reagent"] not in Rgts_Dict):
                Rgts_Dict[rgt["reagent"]]=0
            Rgts_Dict[rgt["reagent"]]+=1

        for rgt in Rgts_Dict.keys():
            if(Rgts_Dict[rgt]>1):
                return "ERROR: Duplicate reagents"

        ########################################
        # Check for duplicate compounds in
        # different compartments, these are 
        # balanced directly.
        #######################################
        Cpds_Coeff_Dict=dict()
        for rgt in Rgts_Array:
            cpd=rgt["compound"]
            if(cpd not in Cpds_Coeff_Dict):
                Cpds_Coeff_Dict[cpd]=0

            #Use float() because you can get real coefficients
            Cpds_Coeff_Dict[cpd]+=float(rgt["coefficient"])

        #Build dict of compounds
        Cpds_Dict=dict()
        for rgt in Rgts_Array:
            rgt["coefficient"]=Cpds_Coeff_Dict[rgt["compound"]]
            Cpds_Dict[rgt["compound"]]=rgt

        ########################################
        # Check for duplicate elements, across
        # all compounds, these are balanced 
        # directly.
        #######################################
        Rxn_Net_Charge=0.0
        Rxn_Net_Mass=dict()
        for cpd in Cpds_Dict.keys():
            cpdAtoms = self.CompoundsHelper.parseFormula(Cpds_Dict[cpd]["formula"])
            if(len(cpdAtoms.keys())==0):
                return "CPDFORMERROR"

            Cpd_Coeff_Charge = float(Cpds_Dict[cpd]["charge"]) * float(Cpds_Dict[cpd]["coefficient"])
            Rxn_Net_Charge += Cpd_Coeff_Charge

            for atom in cpdAtoms.keys():
                Atom_Coeff_Mass = float(cpdAtoms[atom]) * float(Cpds_Dict[cpd]["coefficient"])

                if(atom not in Rxn_Net_Mass.keys()):
                    Rxn_Net_Mass[atom]=0.0

                Rxn_Net_Mass[atom] += Atom_Coeff_Mass

        #Round out tiny numbers that occur because we add/substract floats
        #Threshold of 1e-6 found to capture all these instances without
        #removing actual small differences in mass.
        for atom in Rxn_Net_Mass.keys():
            if(Rxn_Net_Mass[atom] > -1e-6 and Rxn_Net_Mass[atom] < 1e-6):
                Rxn_Net_Mass[atom]=0

        if(Rxn_Net_Charge > -1e-6 and Rxn_Net_Charge < 1e-6):
            Rxn_Net_Charge=0
            
        #Report any imbalance
        Imbalanced_Atoms_Array=list()
        for atom in sorted(Rxn_Net_Mass.keys()):
            if(Rxn_Net_Mass[atom]==0):
                continue

            #Correct for redundant ".0" in floats
            if(str(Rxn_Net_Mass[atom])[-2:] == ".0"):
                Rxn_Net_Mass[atom]=int(round(Rxn_Net_Mass[atom]))

            Imbalanced_Atoms_Array.append(atom+":"+str(Rxn_Net_Mass[atom]))
            
        #Correct for redundant ".0" in floats
        if(str(Rxn_Net_Charge)[-2:] == ".0"):
            Rxn_Net_Charge=int(Rxn_Net_Charge)

        Status=""

        if(len(Imbalanced_Atoms_Array)>0):
            Status = "MI:"+"/".join(Imbalanced_Atoms_Array)

        if(Rxn_Net_Charge != 0):
            if(len(Status)>0):
                Status = "CI:"+str(Rxn_Net_Charge)
            else:
                Status += "|CI:"+str(Rxn_Net_Charge)

        if(Status == ""):
            Status = "OK"

        return Status

    def adjustCompound(self,Rxn_Cpds_Array,Compound,Adjustment,Compartment=0):
        if(Adjustment==0):
            return Rxn_Cpds_Array

        ######################################################################
        # We will always assume to adjust a compound automatically
        # in the compartment indexed as zero, unless otherwise specified.
        # This answers the question of how to handle transporters.
        ######################################################################

        #Check to see if it already exists
        Cpd_Exists=0
        for rgt in Rxn_Cpds_Array:
            if(rgt["compound"]==Compound and rgt["compartment"]==Compartment):
                rgt["coefficient"] += Adjustment
                Cpd_Exists=1

        if(Cpd_Exists != 1):
            rgt_id = Compound+"_"+str(Compartment)+"0"

            Rxn_Cpds_Array.append({"reagent":rgt_id,"coefficient":Adjustment,
                                   "compound":Compound,"compartment":Compartment,"index":0,
                                   "name":self.Compounds_Dict[Compound]["name"],
                                   "formula":self.Compounds_Dict[Compound]["formula"],
                                   "charge":self.Compounds_Dict[Compound]["charge"]})            
        return

    def rebuildReaction(self,Reaction_Dict,Stoichiometry):
        #Assign stoich
        Reaction_Dict["stoichiometry"]=Stoichiometry

        #Build list of "reagents" and "products"
        Rxn_Cpds_Array=self.parseStoich(Stoichiometry)
        Reagents_Array=list()
        Products_Array=list()
        CompoundIDs_Dict=dict()
        for rgt in Rxn_Cpds_Array:
            CompoundIDs_Dict[rgt["compound"]]=1
            if(rgt["coefficient"]>0):
                Products_Array.append(rgt)
            else:
                Reagents_Array.append(rgt)

        RgtsStr_Array=list()
        for rgt in Reagents_Array:
            ID_String = "("+str(abs(rgt["coefficient"]))+") "+rgt["compound"]+"["+str(rgt["index"])+"]"
            RgtsStr_Array.append(ID_String)

        Equation_Array=list()
        Code_Array=list()
        Definition_Array=list()

        Equation_Array.append(" + ".join(RgtsStr_Array))
        Definition_Array.append(" + ".join(RgtsStr_Array))
        Code_Array.append(" + ".join(x for x in RgtsStr_Array if "cpd00067" not in x))
                    
        Code_Array.append("<=>")
        if(Reaction_Dict["direction"] == "="):
            Equation_Array.append("<=>")
            Definition_Array.append("<=>")
        elif(Reaction_Dict["direction"] == "<"):
            Equation_Array.append("<=")
            Definition_Array.append("<=")
        else:
            Equation_Array.append("=>")
            Definition_Array.append("=>")

        PdtsStr_Array=list()
        for rgt in Products_Array:
            ID_String = "("+str(abs(rgt["coefficient"]))+") "+rgt["compound"]+"["+str(rgt["index"])+"]"
            PdtsStr_Array.append(ID_String)

        Equation_Array.append(" + ".join(PdtsStr_Array))
        Definition_Array.append(" + ".join(PdtsStr_Array))
        Code_Array.append(" + ".join(x for x in PdtsStr_Array if "cpd00067" not in x))

        Reaction_Dict["code"]=" ".join(Code_Array)
        Reaction_Dict["equation"]=" ".join(Equation_Array)
        Reaction_Dict["definition"]=" ".join(Definition_Array)
        Reaction_Dict["compound_ids"]=";".join(sorted(CompoundIDs_Dict.keys()))
        
        #Replace ids with names in Definition
        for cpd_id in CompoundIDs_Dict.keys():
            if(cpd_id in Reaction_Dict["definition"]):
                Reaction_Dict["definition"] = Reaction_Dict["definition"].replace(cpd_id,self.Compounds_Dict[cpd_id]["name"])

        return

    def saveReactions(self,Reactions_Dict):
        RxnsRoot = os.path.splitext(self.RxnsFile)[0]
        
        #Print to TSV
        RxnsFile = open(RxnsRoot+".tsv",'w')
        RxnsFile.write("\t".join(self.Headers)+"\n")
        for rxn in sorted(Reactions_Dict.keys()):
            RxnsFile.write("\t".join(str(Reactions_Dict[rxn][header]) for header in self.Headers)+"\n")
        RxnsFile.close()

        #Print to JSON
        RxnsFile = open(RxnsRoot+".json",'w')
        RxnsFile.write(json.dumps(Reactions_Dict, indent=4, sort_keys=True))
        RxnsFile.close()
