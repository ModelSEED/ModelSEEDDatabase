"""Validates BiomassCompounds and Reactions files in Templates folder"""

from csv import DictReader
import os
import argparse


def get_id_set(tsv_path):
    with open(tsv_path) as infile:
        return set([line.split('\t')[0] for line in infile])


def validate_biomass_compounds(path, comp_set):
    with open(path) as infile:
        template_compounds = set()
        for line in DictReader(infile, dialect='excel-tab'):
            template_compounds.add(line['id'])
            if line['linked_compounds'] == 'null':
                continue
            template_compounds.update([x.split(':')[0] for x
                                       in line['linked_compounds'].split('|')])
        return template_compounds - comp_set


def validate_reaction_list(path, rxn_set):
    with open(path) as infile:
        template_rxns = set()
        for line in infile:
            template_rxns.add(line.split('\t')[0])
        return template_rxns - rxn_set

if __name__ == '__main__':
    script_dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser(
        description='Validates BiomassCompounds and Reactions files in Templates folder')
    parser.add_argument('-c', dest='comp_tsv',
                        default=script_dir+'/../Biochemistry/compounds.master.tsv',
                        help='Path to the compounds file')
    parser.add_argument('-r', dest='rxn_tsv',
                        default=script_dir+'/../Biochemistry/reactions.master.tsv',
                        help='Path to the reaction file')
    parser.add_argument('-t', dest='template_dir',
                        default=script_dir+'/../Templates',
                        help='Path to the templates directory')
    args = parser.parse_args()
    comp_ids = get_id_set(args.comp_tsv)
    rxn_ids = get_id_set(args.rxn_tsv)
    exit_code = 0
    for template in os.listdir(args.template_dir):
        if not os.path.isdir(os.path.join(args.template_dir, template)):
            continue
        print("Validating %s template" % template)
        undefined_comps = validate_biomass_compounds('%s/%s/BiomassCompounds.tsv'
                                                   % (args.template_dir, template), comp_ids)
        if undefined_comps:
            print("Undefined Compounds: " + ", ".join(undefined_comps))
            exit_code = 1
        undefined_rxns = validate_reaction_list('%s/%s/Reactions.tsv'
                                                % (args.template_dir, template), rxn_ids)
        if undefined_rxns:
            print("Undefined Reactions: " + ", ".join(undefined_rxns))
            exit_code = 1
    exit(exit_code)
