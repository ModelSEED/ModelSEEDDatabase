import os, re, sys

sys.path.append('../../Libs/Python')
from BiochemPy import Compounds

InChI_Layers = ('c', 'h', 'p', 'q', 'b', 't', 'm', 's')

def parse(inchi, merge_formula=False):
    """
    @param inchi: InChI string
    @param merge_formula: bool, use (not yet implemented) "merge_formulas"
    @return:
    formula string and dictionary of layers where key is layer code and value
    is layer contents
    """
    layer_dict = dict([(x, "") for x in InChI_Layers])

    # special case for proton
    m = re.match('^InChI=1S/p([-+]\d*)', inchi)
    if m:
        layer_dict['p'] = m.group(1)
        return "", layer_dict

    layers = inchi.split("/")[1:]
    formula = layers.pop(0)
    if merge_formula:
        formula = Compounds.mergeFormula(formula)

    for l in layers:
        layer_dict[l[0]] = l[1:]

    return formula, layer_dict

def build(formula, layers, remove=(), merge_formula=False):
    """
    I use 'remove' to strip p, q, and stereochemical layers depending on how I
    want to compare InChI strings
    @param formula: Formula string
    @param layers: layers dictionary
    @param remove: a dictionary of layer codes that have to be removed from InChI string
    @param merge_formula: bool, use (not yet implemented) "merge_formulas"
    @return: InChI string
    """
    if merge_formula:
        formula = Compounds.mergeFormula(formula)
    inchi = "/".join(["InChI=1S"]+[formula]+[layers[x] for x in InChI_Layers
                                             if layers[x] and x not in remove])
    # if no valid layers return blank string
    return inchi if len(inchi) > 8 else ""


def charge(q_layer,p_layer):
    """
    @param q_layer: q layer string
    @param p_layer: p layer string
    @return: charge
    """
    global_charge = 0
    if(q_layer != ""):
        q_components = q_layer.split(';')
        for q_component in q_components:
            if(q_component!=""):
                (multiplier,charge) = (1,0)

                #Match for multiplier
                m = re.match('^(\d+)\*(.+)$', q_component)
                if m:
                    multiplier = int(m.group(1))
                    charge = int(m.group(2))
                else:
                    charge = int(q_component)

                global_charge += ( multiplier * charge )

    #Proton layers have never had multiple components
    #So this is an explicit warning, using it directly 
    #will raise an error
    if(";" in p_layer):
        print("Warning: multiple components in mobile proton layer")

    #protons have positive charge
    if(p_layer != ''):
        global_charge += int(p_layer)

    return global_charge

def adjust_protons(formula, protons):
    """
    @param formula: chemical formula as string
    @param protons: number of hydrogens to add/remove as intager
    @return: new formula as string
    """
    if not protons:
        return (formula,"")
    protons = int(protons)
    Notes = ""
    #The whole function assumes that there is a single formula string
    #If the formula can be broken into components, it must first be merged
    #This is because the proton layer only ever has a single component
    if(len(formula.split('.'))>1):
        print("Error: you must merge the formula components into a single formula string")
        print("You can do so using Compounds.mergeFormula()")
        return formula,"Unadjustable due to multiple components"

    atoms = Compounds.parseFormula(formula)
    if "H" in atoms:
        atoms['H'] += protons
        if atoms['H'] < 0:
            Notes = 'Too Many Protons adjusted!'
        if atoms['H'] == 0:
            del atoms['H']
    elif(len(atoms)==0):
        #special case for the proton
        atoms['H']=protons

    formula = Compounds.buildFormula(atoms)
    return (formula, Notes)
