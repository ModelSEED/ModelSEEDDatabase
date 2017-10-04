import re
from Libs.Python.BiochemPy.Compounds import Compounds

InChI_Layers = ('c', 'h', 'p', 'q', 'b', 't', 'm', 's')


def parse(inchi, merge_formula=False):
    """
    @param inchi: InChI string
    @param merge_formula: bool, use (not yet implemented) "merge_formulas"
    @return:
    formula string and dictionary of layers where key is layer code and value
    is layer contents
    """
    # special case for proton
    m = re.match('^InChI=1S/p([-+]\d*)', inchi)
    if m:
        return "", {'p': m.group(1)}
    layers = inchi.split("/")[1:]
    formula = layers.pop(0)
    if merge_formula:
        formula = Compounds.mergeFormula(formula)
    layer_dict = dict([(x, "") for x in InChI_Layers])
    for l in layers:
        layer_dict[l[0]] = l[1:]
    return formula, layers


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


def charge(q_str,p_str):
    """
    Not sure how this one works
    @param q_str: q layer string
    @param p_str: p layer string
    @return: charge
    """
    raise NotImplementedError


def adjust_protons(formula, protons):
    """

    @param formula: chemical formula as string
    @param protons: number of hydrogens to add/remove as intager
    @return: new formula as string
    """
    if not protons:
        return formula
    formula_parts = formula.split('.')
    for i, component in enumerate(formula_parts):
        atoms = Compounds.parseFormula(component)
        if "H" in atoms:
            atoms['H'] += protons
            if atoms['H'] < 0:
                print('ERROR: Too Many Protons adjusted in formula!')
            if atoms['H'] == 0:
                del atoms['H']
            formula_parts[i] = Compounds.hill_sorted(atoms)
            break
    return '.'.join(formula_parts)
