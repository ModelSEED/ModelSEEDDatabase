#!/usr/bin/env python
from BiochemPy import Reactions
from math import log
reactions_helper = Reactions()
reactions_dict = reactions_helper.loadReactions()

#Constants
TEMPERATURE=298.15
GAS_CONSTANT=0.0019858775
RT_CONST=TEMPERATURE*GAS_CONSTANT
FARADAY = 0.023061 # kcal/vol gram divided by 1000?

# max and min values refer to range of intracellular concentrations
(cell_max,cell_min,cell_conc)=(0.02,0.00001,0.001)

#Phosphates
phosphate_ids=("cpd00002", #ATP
               "cpd00008", #ADP
               "cpd00018", #AMP
               "cpd00009", #Pi
               "cpd00012") #PPi

#Low energy compounds
#taken from MFAToolkit/Parameters/Defaults.txt
low_energy_cpds=("cpd00011", #CO2
                 "cpd00013", #NH3
                 "cpd11493", #ACP
                 "cpd00009", #Pi
                 "cpd00012", #Ppi
                 "cpd00010", #CoA
                 "cpd00449", #Dihydrolipoamide
                 "cpd00242") #HCO3

for rxn in sorted(reactions_dict.keys()):
    
    #defaults
    thermoreversibility = "?"
    proton_cpt_dict = dict()

    if(reactions_dict[rxn]['status'] == "EMPTY"):
        print(rxn+":EMPTY")
        reactions_dict[rxn]['reversibility']="?"
        continue

    rxn_dg = reactions_dict[rxn]['deltag']
    rxn_dge = reactions_dict[rxn]['deltagerr']

    if(rxn_dg == 10000000 or rxn_dg is None):
        print(rxn+":Incomplete:"+reactions_dict[rxn]["reversibility"])
        reactions_dict[rxn]['reversibility']="?"
        continue

    #Calculate MdeltaG
    (rct_min,rct_max)=(0.0,0.0)
    (pdt_min,pdt_max)=(0.0,0.0)

    #Calculate mMdeltaG
    rgt_sum=0.0

    #Capture specific compounds for heuristics
    Phosphates=dict()
    for rgt in reactions_dict[rxn]['stoichiometry'].split(';'):
        (coeff,cpd,cpt,idx,name)=rgt.split(":",maxsplit=4)
        coeff=float(coeff)

        if(cpd == 'cpd00067'):
            proton_cpt_dict[cpt]=1

        #Find phosphates
        for cpd in phosphate_ids:
            if(cpd in rgt):
                if(cpd not in Phosphates):
                    Phosphates[cpd]=0.0
                Phosphates[cpd]+=coeff

        #ignore protons and water for following computation
        if(cpd == 'cpd00067' or cpd == 'cpd00001'):
            continue

        #Here we can change accordingly to compartments
        #This section for MdeltaG under concentration range
        (cpt_max,cpt_min)=(cell_max,cell_min)
        if(coeff<0):
            rct_min += (coeff*log(cpt_min))
            rct_max += (coeff*log(cpt_max))
        else:
            pdt_min += (coeff*log(cpt_min))
            pdt_max += (coeff*log(cpt_max))

        #This section for mMdeltaG under fixed concentration
        local_conc=cell_conc
        if(cpd == 'cpd00011'): #CO2
            local_conc=0.0001
        elif(cpd == 'cpd00007' or cpd == 'cpd11640'): #O2 && H2
            local_conc=0.000001;
        rgt_sum += (coeff*log(local_conc))

    #for future reference
    rxn_dg_transport = 0.0
    
    stored_max=rxn_dg+rxn_dg_transport+rxn_dge
    stored_min=rxn_dg+rxn_dg_transport+rxn_dge

    stored_max+=(RT_CONST*pdt_max)+(RT_CONST*rct_min)
    stored_min+=(RT_CONST*pdt_min)+(RT_CONST*rct_max)

    if(stored_max < 0):
        thermoreversibility = ">"
        print(rxn+":MdeltaG:>:{0:.4f}".format(stored_max))
        reactions_dict[rxn]['reversibility']=">"
        continue
    if(stored_min > 0):
        thermoreversibility = "<"
        print(rxn+":MdeltaG:<:{0:.4f}".format(stored_min))
        reactions_dict[rxn]['reversibility']="<"
        continue

    #Do heuristics
    #1: ATP hydrolysis transport
    #1a: ATP Synthase is reversible, but cannot involve any other compound, and can only transport protons
    is_atp_synthase=False
    if(reactions_dict[rxn]['is_transport']==1 and len(proton_cpt_dict.keys())>1):

        cpds_cpts_dict=dict()
        #Collect compound compartments
        for rgt in reactions_dict[rxn]['stoichiometry'].split(';'):
            (coeff,cpd,cpt,idx,name)=rgt.split(":",maxsplit=4)
            coeff=float(coeff)
        
            if(cpd not in cpds_cpts_dict):
                cpds_cpts_dict[cpd]=list()
            cpds_cpts_dict[cpd].append(cpt)

        #defaults
        is_atp_synthase=True
        for cpd in cpds_cpts_dict.keys():
            #Must not contain reactants not in ATP Synthase
            if(cpd != 'cpd00002' and cpd != 'cpd00008' and cpd != 'cpd00009' and cpd != 'cpd00001' and cpd != 'cpd00067'):
                is_atp_synthase = False

        #Must contain _all_ five reactants in ATP Synthase
        if(len(cpds_cpts_dict.keys())!=5):
            is_atp_synthase = False

        #Only protons are transported
        for cpd in cpds_cpts_dict.keys():
            if(len(cpds_cpts_dict[cpd])==2 and cpd != 'cpd00067'):
                is_atp_synthase = False

    if(is_atp_synthase is True):
        thermoreversibility="="
        reactions_dict[rxn]['reversibility']="="
        print(rxn+":ATPS:"+thermoreversibility)
        continue

    #1b: Find ABC Transporters (but not ATP Synthase)
    if(reactions_dict[rxn]['is_transport']==1 and 'cpd00002' in Phosphates):
        if(Phosphates['cpd00002']<0):
            thermoreversibility=">"
        elif(Phosphates['cpd00002']>0):
            thermoreversibility="<"
        else:
            #If zero, then itself ATP is transported
            #I manually reviewed these, these are not chemical reactions
            thermoreversibility="="
        
        reactions_dict[rxn]['reversibility']=thermoreversibility
        print(rxn+":ABCT:"+str(Phosphates['cpd00002'])+":"+thermoreversibility)

    #2: Calculate and evaluate mMdeltaG
    mMdeltaG=rxn_dg+(RT_CONST*rgt_sum);
    if(mMdeltaG >= -2.0 and mMdeltaG <= 2.0):
        thermoreversibility="="
        reactions_dict[rxn]['reversibility']="="
        print(rxn+":mMdeltaG:=:{0:.4f}".format(mMdeltaG))
        continue

    #3: Calculate low energy points
    low_energy_points = 0

    #3a: Find minimum phosphate-related coefficient
    min_coeff = 10000000
    if('cpd00002' in Phosphates and len(Phosphates.keys())>2):
        print(":Phos",Phosphates)
        for pho in Phosphates.keys():
            if(Phosphates[pho]<min_coeff):
                min_coeff = Phosphates[pho]

    if(min_coeff != 10000000):
        low_energy_points-=(abs(min_coeff))
    
    #3b:Find other low energy compounds
    for rgt in reactions_dict[rxn]['stoichiometry'].split(';'):
        (coeff,cpd,cpt,idx,name)=rgt.split(":",maxsplit=4)
        coeff=float(coeff)
        
        if(cpd in low_energy_cpds):
            low_energy_points-=coeff

    #Evaluate low energy
    if((low_energy_points*mMdeltaG) > 2 and mMdeltaG < 0):
        thermoreversibility = ">"
        reactions_dict[rxn]['reversibility']=">"
        print(rxn+":lowE:=:{0:.4f}:".format(mMdeltaG)+str(low_energy_points)+":>")
        continue
    elif((low_energy_points*mMdeltaG) > 2 and mMdeltaG > 0):
        thermoreversibility = "<"
        reactions_dict[rxn]['reversibility']="<"
        print(rxn+":lowE:=:{0:.4f}:".format(mMdeltaG)+str(low_energy_points)+":<")
        continue

    thermoreversibility = "="
    reactions_dict[rxn]['reversibility']="="
    print(rxn+":Defaut:=")
reactions_helper.saveReactions(reactions_dict)
