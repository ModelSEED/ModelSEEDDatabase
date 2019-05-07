from . import Schemas as schemas
from BiochemPy import Reactions
reactions_helper = Reactions()

from .error_reporting import find_new_errors, report_errors
from jsonschema import Draft4Validator
from collections import defaultdict, Counter
import re
import argparse
import json


def parse_rxn(_rxn):
    # returns a tuple of dicts with compounds as keys and stoich as values
    _rxn = _rxn.translate(str.maketrans("", "", "()"))
    return [dict([compound[:-3].split()[::-1] for compound in half.split(' + ')])
            if half else {} for half in re.split(' <?=>? ', _rxn)]


def get_atom_count(compoundDict, complist):
    atom_counts = Counter()
    for id, stoich in complist.items():
        for comp in compoundDict:
            if(id != comp['id']):
                continue
            if comp['formula'] is None or comp['formula'] == 'null':
                continue
            for pair in re.findall('([A-Z][a-z]?)(\d*)', comp['formula']):
                if not pair[1]:
                    atom_counts[pair[0]] += float(stoich)
                else:
                    atom_counts[pair[0]] += float(pair[1]) * float(stoich)
    atom_counts = Counter({k: round(v, 1) for k, v in atom_counts.items()})
    return atom_counts


def validate_schema(_rxns, verbose):
    v = Draft4Validator(schemas.reactions)
    schema_errors = sorted(v.iter_errors(rxns), key=lambda x: x.path)
    if verbose:
        for e in schema_errors:
            print("Schema Validation ERROR at reaction[{}]: {}"
                  .format("][".join(e.absolute_path), e.message))
    return len(schema_errors)


def check_dups(_rxns, verbose, unique_fields=('id', 'stoichiometry')):
    # build a nested dict for uniqueness checking
    # {field: {value: [ids_with_this_value]}}
    unique_values = dict([(key, defaultdict(list)) for key in unique_fields])
    for rxn in _rxns:
        if rxn['is_obsolete']:
            continue
        if 'EMPTY' in rxn['status']:
            continue

        for key in unique_values:
            unique_values[key][rxn[key]].append(rxn['id'])

    # if the unique_values dict for a field has more than one id, it's not unique
    duplicated = defaultdict(int)
    for value_type in unique_values:
        for key, ids in unique_values[value_type].items():
            if len(ids) > 1 and key != 'null':
                if verbose:
                    print("Duplicate {}: {} in {}".format(value_type, key, ids))
                duplicated["duplicated_"+value_type] += len(ids)-1
    return duplicated


def check_compounds(_rxns, verbose, compound_loc='./Biochemistry/compounds.json'):
    err = defaultdict(int)
    compounds = json.load(open(compound_loc))
    obsolete_comps = set(comp['id'] for comp in compounds
                         if comp['is_obsolete'] == '1')
    for rxn in _rxns:
        if rxn['is_obsolete'] == "1":
            continue
        comp_ids = set(rxn['compound_ids'].split(';'))
        if comp_ids & obsolete_comps:
            if verbose:
                print('Obsolete compounds in {}: {}'
                      .format(rxn['id'], comp_ids & obsolete_comps))
            err['obsolete_comps'] += 1
        try:
            reactants, products = parse_rxn(rxn['equation'])
            rxn_cpds_array=reactions_helper.parseStoich(rxn["stoichiometry"])
        except ValueError as e:
            print("Unable to parse {}: {}".format(rxn['equation'], e))
            err['invalid_equation'] += 1

        try:
            new_status = reactions_helper.balanceReaction(rxn_cpds_array)
#            reactant_atoms = get_atom_count(compounds, reactants)
#            product_atoms = get_atom_count(compounds, products)
        except KeyError as e:
            print('Invalid id {} in equation {}'.format(rxn['id'], e))
            err['invalid_equation'] += 1

#        if reactant_atoms - product_atoms or product_atoms - reactant_atoms:
        if 'OK' not in new_status:
            err['unbalanced_reactions'] += 1
            if rxn['status'] == 'OK':
                if verbose:
                    print('Unbalanced reaction marked OK {}\n{}\n{}'
                          .format(rxn['id'], reactant_atoms, product_atoms))
                err['unbalanced_marked_OK'] += 1
        if rxn['status'] == "EMPTY":
            print("Empty Reaction: " + rxn['id'])
            continue
        if comp_ids ^ set(list(reactants.keys()) + list(products.keys())):
            if verbose:
                print('"compound_ids" and "equation" are inconsistant in {}:\n{}\n{}'
                      .format(rxn['id'], rxn['compound_ids'], rxn['code']))
            err['inconsistent_equation_compids'] += 1
    return err


if __name__ == "__main__":
    # Parse options.
    parser = argparse.ArgumentParser()
    parser.add_argument('rxnfile', help='path to reactions file', action='store')
    parser.add_argument('-d', '--print-duplicates', help='Print all duplicated values', action='store_true', dest='print_dups', default=False)
    parser.add_argument('-s', '--print-schema', help='Print all schema violations', action='store_true', dest='print_schema', default=False)
    parser.add_argument('-c', '--print-compound', help='Print all compound related violations', action='store_true', dest='print_compound', default=False)
    args = parser.parse_args()
    rxns = json.load(open(args.rxnfile))
    errors = dict()
    errors['schema_violations'] = validate_schema(rxns, args.print_schema)
    errors.update(check_dups(rxns, args.print_dups))
    errors.update(check_compounds(rxns, args.print_compound))
    print("\nError counts")
    for error_type, count in errors.items():
        print("\t{}: {}".format(error_type, count))
    new_errors = find_new_errors('reactions', errors)
    report_errors('reactions', errors)
    print("New errors detected: " + ", ".join(new_errors))
    exit(len(new_errors))
