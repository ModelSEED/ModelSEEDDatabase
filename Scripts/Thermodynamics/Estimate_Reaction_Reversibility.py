#!/usr/bin/env python
import sys
from math import log

# Constants
TEMPERATURE = 298.15
GAS_CONSTANT = 0.0019858775
RT_CONST = TEMPERATURE * GAS_CONSTANT
FARADAY = 0.023061  # kcal/vol gram divided by 1000?

# max and min values refer to range of intracellular concentrations
(cell_max, cell_min, cell_conc) = (0.02, 0.00001, 0.001)

# Sentinel value marking a missing energy estimate
SENTINEL = 10000000

# Phosphates
phosphate_ids = ("cpd00002",  # ATP
                 "cpd00008",  # ADP
                 "cpd00018",  # AMP
                 "cpd00009",  # Pi
                 "cpd00012")  # PPi

# Low energy compounds
# taken from MFAToolkit/Parameters/Defaults.txt
low_energy_cpds = ("cpd00011",  # CO2
                   "cpd00013",  # NH3
                   "cpd11493",  # ACP
                   "cpd00009",  # Pi
                   "cpd00012",  # Ppi
                   "cpd00010",  # CoA
                   "cpd00449",  # Dihydrolipoamide
                   "cpd00242")  # HCO3


def _estimate_core(rxn_obj, rxn_dg, rxn_dge):
    """Core thermodynamic-direction heuristic for a single deltaG estimate.

    Given a reaction object and one usable (deltaG, deltaG-error) estimate,
    return an (operator, reason) tuple where operator is one of '>', '<', '='.
    The reason string mirrors the first column written to the reversibility
    report. Callers must screen out EMPTY / missing / sentinel estimates before
    calling this (see reversibility_from_energy).

    The logic is preserved verbatim from the original Estimate_Reaction_
    Reversibility.py so that the value returned for the estimate currently
    driving a reaction's deltag matches that reaction's canonical reversibility.
    """

    # Calculate MdeltaG
    (rct_min, rct_max) = (0.0, 0.0)
    (pdt_min, pdt_max) = (0.0, 0.0)

    # Calculate mMdeltaG
    rgt_sum = 0.0

    # Capture specific compounds for heuristics
    proton_cpt_dict = dict()
    phosphates = dict()
    for rgt in rxn_obj['stoichiometry']:
        cpd = rgt['compound']
        cpt = rgt['compartment']
        coeff = float(rgt['coefficient'])

        if(cpd == 'cpd00067'):
            proton_cpt_dict[cpt] = 1

        # Find phosphates
        for cpd in phosphate_ids:
            if(cpd in rgt):
                if(cpd not in phosphates):
                    phosphates[cpd] = 0.0
                phosphates[cpd] += coeff

        # ignore protons and water for following computation
        if(cpd == 'cpd00067' or cpd == 'cpd00001'):
            continue

        # Here we can change accordingly to compartments
        # This section for MdeltaG under concentration range
        (cpt_max, cpt_min) = (cell_max, cell_min)
        if(coeff < 0):
            rct_min += (coeff * log(cpt_min))
            rct_max += (coeff * log(cpt_max))
        else:
            pdt_min += (coeff * log(cpt_min))
            pdt_max += (coeff * log(cpt_max))

        # This section for mMdeltaG under fixed concentration
        local_conc = cell_conc
        if(cpd == 'cpd00011'):  # CO2
            local_conc = 0.0001
        elif(cpd == 'cpd00007' or cpd == 'cpd11640'):  # O2 && H2
            local_conc = 0.000001
        rgt_sum += (coeff * log(local_conc))

    # for future reference
    rxn_dg_transport = 0.0

    stored_max = rxn_dg + rxn_dg_transport + rxn_dge
    stored_min = rxn_dg + rxn_dg_transport - rxn_dge

    stored_max += (RT_CONST * pdt_max) + (RT_CONST * rct_min)
    stored_min += (RT_CONST * pdt_min) + (RT_CONST * rct_max)

    if(stored_max < 0):
        return (">", "MdeltaG(Max): {0:.2f}".format(stored_max))

    if(stored_min > 0):
        return ("<", "MdeltaG(Min): {0:.2f}".format(stored_min))

    # Do heuristics
    # 1: ATP hydrolysis transport
    # 1a: ATP Synthase is reversible, but cannot involve any other compound, and can only transport protons
    is_atp_synthase = False
    if(rxn_obj['is_transport'] == 1 and len(proton_cpt_dict.keys()) > 1):
        cpds_cpts_dict = dict()
        # Collect compound compartments
        for rgt in rxn_obj['stoichiometry']:
            cpd = rgt['compound']
            cpt = rgt['compartment']
            coeff = float(rgt['coefficient'])

            if(cpd not in cpds_cpts_dict):
                cpds_cpts_dict[cpd] = list()
            cpds_cpts_dict[cpd].append(cpt)

        # defaults
        is_atp_synthase = True
        for cpd in cpds_cpts_dict.keys():
            # Must not contain reactants not in ATP Synthase
            if(cpd != 'cpd00002' and cpd != 'cpd00008' and cpd != 'cpd00009' and cpd != 'cpd00001' and cpd != 'cpd00067'):
                is_atp_synthase = False

        # Must contain _all_ five reactants in ATP Synthase
        if(len(cpds_cpts_dict.keys()) != 5):
            is_atp_synthase = False

        # Only protons are transported
        for cpd in cpds_cpts_dict.keys():
            if(len(cpds_cpts_dict[cpd]) == 2 and cpd != 'cpd00067'):
                is_atp_synthase = False

    if(is_atp_synthase is True):
        return ("=", "ATPS")

    # 1b: Find ABC Transporters (but not ATP Synthase)
    if(rxn_obj['is_transport'] == 1 and 'cpd00002' in phosphates):

        thermoreversibility = "="

        if(phosphates['cpd00002'] < 0):
            thermoreversibility = ">"
        elif(phosphates['cpd00002'] > 0):
            thermoreversibility = "<"
        else:
            # If zero, then itself ATP is transported
            # I manually reviewed these, these are not chemical reactions
            pass

        return (thermoreversibility, "ABCT: " + str(phosphates['cpd00002']))

    # 2: Calculate and evaluate mMdeltaG
    mMdeltaG = rxn_dg + (RT_CONST * rgt_sum)
    if(mMdeltaG >= -2.0 and mMdeltaG <= 2.0):
        return ("=", "mMdeltaG: {0:.2f}".format(mMdeltaG))

    # 3: Calculate low energy points
    low_energy_points = 0

    # 3a: Find minimum phosphate-related coefficient
    min_coeff = 10000000
    if('cpd00002' in phosphates and len(phosphates.keys()) > 2):
        for pho in phosphates.keys():
            if(phosphates[pho] < min_coeff):
                min_coeff = phosphates[pho]

    if(min_coeff != 10000000):
        low_energy_points -= (abs(min_coeff))

    # 3b:Find other low energy compounds
    for rgt in rxn_obj['stoichiometry']:
        cpd = rgt['compound']
        cpt = rgt['compartment']
        coeff = float(rgt['coefficient'])

        if(cpd in low_energy_cpds):
            low_energy_points -= coeff

    # Evaluate low energy
    if((low_energy_points * mMdeltaG) > 2 and mMdeltaG < 0):
        return (">", "lowE: {0:.2f}".format(mMdeltaG) + ":" + str(low_energy_points))

    elif((low_energy_points * mMdeltaG) > 2 and mMdeltaG > 0):
        return ("<", "lowE: {0:.2f}".format(mMdeltaG) + ":" + str(low_energy_points))

    return ("=", "default")


def reversibility_from_energy(rxn_obj, rxn_dg, rxn_dge):
    """Return the thermodynamic-direction operator implied by a single deltaG
    estimate for a reaction.

    This is the per-estimate entry point shared by the canonical reversibility
    script and the per-method thermodynamics records (Group contribution,
    eQuilibrator, dGPredictor). It returns '?' when there is no usable estimate
    (EMPTY reaction, missing value, or the 10000000 sentinel) and otherwise one
    of '>', '<', '='.
    """

    if(isinstance(rxn_obj, dict) and rxn_obj.get('status') == "EMPTY"):
        return "?"

    # deltag may be stored as a numeric string; coerce as the canonical
    # reversibility script does. Anything non-coercible, None, or the 10000000
    # sentinel is treated as "no usable estimate".
    if(rxn_dg is None or isinstance(rxn_dg, bool)):
        return "?"
    try:
        rxn_dg = float(rxn_dg)
    except (TypeError, ValueError):
        return "?"
    if(rxn_dg == SENTINEL):
        return "?"

    try:
        rxn_dge = 0.0 if (rxn_dge is None or isinstance(rxn_dge, bool)) else float(rxn_dge)
    except (TypeError, ValueError):
        rxn_dge = 0.0

    (operator, _reason) = _estimate_core(rxn_obj, rxn_dg, rxn_dge)
    return operator


if __name__ == '__main__':

    sys.path.append('../../Libs/Python/')
    from BiochemPy import Reactions
    reactions_helper = Reactions()
    reactions_dict = reactions_helper.loadReactions()

    DB_Level = ''
    if(len(sys.argv) > 1 and (sys.argv[1] == 'EQ' or sys.argv[1] == 'GC')):
        DB_Level = sys.argv[1]

    reversibility_report = dict()
    for rxn in sorted(reactions_dict.keys()):

        rxn_obj = reactions_dict[rxn]

        # defaults
        thermoreversibility = "?"

        if(rxn_obj['status'] == "EMPTY"):

            thermoreversibility = "?"
            reversibility_report[rxn] = ["Empty", rxn_obj["reversibility"], thermoreversibility]
            rxn_obj['reversibility'] = thermoreversibility

            continue

        rxn_dg = rxn_obj['deltag']
        rxn_dge = rxn_obj['deltagerr']

        # Here, if I'm specifying either GC or EQ,
        # Then I want to check that I should estimate for this reaction
        # (I.e. either "GCC" or "EQC")
        # Otherwise its labeled as incomplete
        DB_Rxn = True
        if(len(DB_Level) > 0):
            DB_Rxn = False
            for entry in rxn_obj["notes"]:
                if(DB_Level in entry and (entry == "GCC" or entry == "EQU")):
                    DB_Rxn = True

        # Preserve dev's coercion fix: deltag may be a numeric string.
        if(rxn_dg is not None):
            rxn_dg = float(rxn_dg)

        if(rxn_dg == SENTINEL or rxn_dg is None or DB_Rxn is False):

            thermoreversibility = "?"
            status = "Incomplete"

            # Here, if using EQ, but incomplete/not-updated reaction
            # Can still fall back onto GC if complete by GC standards

            if(DB_Level == "EQ" and "GCC" in rxn_obj['notes']):
                thermoreversibility = rxn_obj["reversibility"]
                status += " (GCC)"

            reversibility_report[rxn] = [status, rxn_obj["reversibility"], thermoreversibility]
            rxn_obj['reversibility'] = thermoreversibility

            continue

        (thermoreversibility, reason) = _estimate_core(rxn_obj, rxn_dg, rxn_dge)
        reversibility_report[rxn] = [reason, rxn_obj["reversibility"], thermoreversibility]
        rxn_obj['reversibility'] = thermoreversibility

    file_name = "Estimated_Reaction_Reversibility_Report"
    if(len(DB_Level) > 0):
        file_name += "_" + DB_Level
    file_name += ".txt"
    with open(file_name, "w") as fh:
        for rxn in sorted(reversibility_report):
            report_array = list(reversibility_report[rxn])
            if(DB_Level == "GC"):
                del(report_array[1])
            fh.write(rxn + "\t" + "\t".join(report_array) + "\n")
    fh.close()

    print("Saving reactions")
    reactions_helper.saveReactions(reactions_dict)
