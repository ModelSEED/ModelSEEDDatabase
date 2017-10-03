import argparse
from ast import literal_eval
from csv import DictReader
import json


def parse_complexes(fp):
    complexes = []
    r = DictReader(open(fp), dialect='excel-tab')
    for line in r:
        line['complexroles'] = []
        line['confidence'] = float(line['confidence'])
        for role in line['roles'].split("|"):
            if role == "null":
                continue
            sp_role = role.split(';')
            line['complexroles'].append({
                'templaterole_ref': '~/roles/id/' + sp_role[0],
                'optional_role': int(sp_role[2]),
                'triggering': 1
            })
        del line['roles']
        complexes.append(line)
    return complexes


def parse_roles(fp):
    roles = []
    r = DictReader(open(fp), dialect='excel-tab')
    for line in r:
        line['features'] = line['features'].split(';')
        line['aliases'] = line['aliases'].split(';')
        roles.append(line)
    return roles


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('template', help='path to template file to edit',
                        action='store')
    parser.add_argument('-a', '--update_all', help='update_all_data_types',
                        action='store_true', default=False)
    parser.add_argument('-R', '--update_roles',
                        help='update roles with Annotations/Roles.tsv',
                        action='store_true', default=False)
    parser.add_argument('-C', '--update_complexes',
                        help='update complexes with Annotations/Roles.tsv',
                        action='store_true', default=False)
    parser.add_argument('-r', '--update_reactions',
                        help='',
                        action='store_true', default=False)

    parser.add_argument('-o', '--other', help='provide a string parseable as '
                        'dict to update other properties', default="")
    args = parser.parse_args()
    if args.update_all:
        args.update_complexes = True
        args.update_roles = True
        args.update_reactions = True
    template_obj = json.load(open(args.template))
    if args.other:
        other_dict = literal_eval(args.other)
        if not isinstance(other_dict, dict) or any(
                [x not in template_obj for x in other_dict.keys()]):
            raise ValueError('Unable to parse other dict')
        template_obj.update(other_dict)
    if args.update_roles:
        template_obj['roles'] = parse_complexes('Annotations/Roles.tsv')
    if args.update_complexes:
        template_obj['complexes'] = parse_complexes('Annotations/Complexes.tsv')
    out_fp = 'Objects/{}.json'.format(template_obj['name'].split('/')[-1])
    json.dump(template_obj, open(out_fp, 'w'))