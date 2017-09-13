"""Validates BiomassCompounds and Reactions files in Templates folder"""

from csv import DictReader
import os
import sys
import argparse

imbalenced = {}

def get_id_set(tsv_path):
    with open(tsv_path) as infile:
        global imbalenced
        ids = set()
        obs = {}
        for line in DictReader(infile, dialect='excel-tab'):
            if 'is_obsolete' in line and int(line['is_obsolete']):
                obs[line['id']] = line.get('linked_reaction', "").split(';')[0]
                continue
            if 'status' in line and ("MI" in line['status']):
                imbalenced[line['id']] = ", ".join([line['id'], line['name'], line['code'], line['status']])
                continue
            ids.add(line['id'])
    return ids, obs


def get_id_dict(tsv_path):
    with open(tsv_path) as infile:
        id_dict = {}
        for line in DictReader(infile, dialect='excel-tab'):
            id_dict[line['id']] = line
    return id_dict


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


def validate_reaction_list(path, rxn_set, complex_set):
    with open(path) as infile:
        template_rxns = set()
        template_complexes = set()
        for line in DictReader(infile, dialect='excel-tab'):
            template_rxns.add(line['id'])
            template_complexes.update(line['complexes'].strip().split("|"))
        return template_rxns - rxn_set, template_complexes - complex_set


def remove_ids(path, ids):
    ids = set(ids)
    txt = open(path).readlines()
    with open(path, 'w') as outfile:
        outfile.writelines([x for x in txt if (x.split('\t')[0] not in ids)])

def update_obsolete(path, obs_rxns):
    txt = open(path).readlines()
    with open(path, 'w') as outfile:
        for line in txt:
            rid = line.split('\t')[0]
            if rid in obs_rxns and obs_rxns[rid]:
                outfile.write(line.replace(rid, obs_rxns[rid]))
            else:
                outfile.write(line)


if __name__ == '__main__':
    script_dir = os.path.dirname(__file__)
    parser = argparse.ArgumentParser(
        description='Validates BiomassCompounds and Reactions files in Templates folder')
    parser.add_argument('-c', dest='comp_tsv',
                        default=script_dir+'/../../Biochemistry/compounds.tsv',
                        help='Path to the compounds file')
    parser.add_argument('-r', dest='rxn_tsv',
                        default=script_dir+'/../../Biochemistry/reactions.tsv',
                        help='Path to the reaction file')
    parser.add_argument('-C', dest='complex_tsv',
                        default=script_dir+'/../../Annotations/Complexes.tsv',
                        help='Path to the complexes file')
    parser.add_argument('-t', dest='template_dir',
                        default=script_dir+'/../../Templates',
                        help='Path to the templates directory')
    parser.add_argument('-d', dest='delete',
                        default=False, action='store_true',
                        help='Delete reactions which are invalid')
    parser.add_argument('-u', dest='update',
                        default=False, action='store_true',
                        help='Update reactions which are invalid')
    args = parser.parse_args()
    comp_ids, obs_comps = get_id_set(args.comp_tsv)
    rxn_ids, obs_rxns = get_id_set(args.rxn_tsv)
    complex_ids = get_id_set(args.complex_tsv)[0] | {'universal', 'null'}
    exit_code = 0
    #The following reactions are inherently unbalanced so we give them a pass
    whitelist_rxns = {'rxn05296', 'rxn05294', 'rxn11921', 'rxn11922'}
    for template in os.listdir(args.template_dir):
        # Skip validation on Mayo clinic templates
        if template in {'Human', 'Microbial', 'Mycobacteria'}:
            continue
        if not os.path.isdir(os.path.join(args.template_dir, template)):
            continue
        print("Validating %s template" % template)
        undef_comps = validate_biomass_compounds(
            '%s/%s/BiomassCompounds.tsv' % (args.template_dir, template),
            comp_ids)
        undef_rxns, undef_complex = validate_reaction_list(
            '%s/%s/Reactions.tsv' % (args.template_dir, template), rxn_ids,
            complex_ids)

        if undef_comps:
            print("ERROR-%s Invalid Compounds: %s"
                  % (len(undef_comps), ", ".join(undef_comps)))
            exit_code = 1
        undef_rxns -= whitelist_rxns
        if undef_rxns:
            print("ERROR-%s Invalid Reactions: %s"
                  % (len(undef_rxns), ", ".join(undef_rxns)))
            for x in undef_rxns:
                if x in imbalenced:
                    print(imbalenced[x])
            exit_code = 1
        if undef_complex:
            print("ERROR-%s Invalid Complexes: %s"
                  % (len(undef_complex), ", ".join(undef_complex)))
            exit_code = 1
        if args.update and undef_rxns:
            update_obsolete('%s/%s/Reactions.tsv' % (args.template_dir, template), obs_rxns)
        if args.delete and undef_comps:
            remove_ids('%s/%s/BiomassCompounds.tsv' % (args.template_dir, template), undef_comps)
        if args.delete and undef_rxns:
            remove_ids('%s/%s/Reactions.tsv' % (args.template_dir, template), undef_rxns)

    exit(exit_code)
