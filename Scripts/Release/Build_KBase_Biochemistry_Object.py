#!/usr/bin/env python
from collections import defaultdict
import csv
from datetime import datetime
import json
import subprocess

local_head = subprocess.check_output(
            ['git', 'log', '-n', '1', '--pretty=format:%H']).decode("utf-8")


def parse_compounds(fp):
    trans_dict = {
        'abbreviation': ('abbreviation', str),
        'deltag': ('deltaG', float),
        'deltagerr': ('deltaGErr', float),
        'formula': ('formula', str),
        'id': ('id', str),
        'name': ('name', str),
        'mass': ('mass', float),
        'charge': ('defaultCharge', int),
        'is_cofactor': ('isCofactor', int),
    }
    source = json.load(open(fp))
    compounds = []
    for in_comp in source.values():
        if in_comp['is_obsolete']:
            continue
        out_comp = {"unchargedFormula": "", 'cues': {}}
        for in_id, out_tup in trans_dict.items():
            try:
                out_comp[out_tup[0]] = out_tup[1](in_comp[in_id])
            except ValueError:
                continue
        out_comp['pkas'] = dict([(int(x.split(':')[0]), [float(x.split(':')[1])])
                                 for x in in_comp['pka'].split(';') if x])
        out_comp['pkbs'] = dict([(int(x.split(':')[0]), [float(x.split(':')[1])])
                                 for x in in_comp['pkb'].split(';') if x])
        compounds.append(out_comp)
    return compounds


def parse_reactions(fp):
    # mapping of field in old to field in new and type
    trans_dict = {
        'abbreviation': ('abbreviation', str),
        'deltag': ('deltaG', float),
        'deltagerr': ('deltaGErr', float),
        'direction': ('direction', str),
        'id': ('id', str),
        'name': ('name', str),
        'reversibility': ('thermoReversibility', str),
        'status': ('status', str),
    }
    source = json.load(open(fp))
    reactions = []
    for in_rxn in source.values():
        if in_rxn['is_obsolete']:
            continue
        out_rxn = {"defaultProtons": 0, 'cues': {}}
        for in_id, out_tup in trans_dict.items():
            try:
                out_rxn[out_tup[0]] = out_tup[1](in_rxn[in_id])
            except ValueError:
                continue
        out_rxn['reagents'] = []
        for comp in in_rxn['stoichiometry'].split(';'):
            sp_comp = comp.split(":")
            try:
                out_rxn['reagents'].append({
                    "coefficient": float(sp_comp[0]),
                    "compartment_ref": "~/compartments/id/" +
                                       ["c", "e", 'p'][int(sp_comp[2])],
                    "isCofactor": 0,
                    "compound_ref": "~/compounds/id/" + sp_comp[1]
                })
            except (IndexError, ValueError):
                in_rxn['stoichiometry']
        reactions.append(out_rxn)
    return reactions


def parse_aliases(fp):
    r = csv.DictReader(open(fp), dialect='excel-tab')
    if r.fieldnames != ['MS ID', 'Old MS ID', 'External ID', 'Source']:
        raise ValueError("Invalid aliases table format")
    aliases = defaultdict(lambda: defaultdict(list))
    for line in r:
        aliases[line['MS ID']][line['Source']].append(line['External ID'])
    return aliases


if __name__ == "__main__":
    biochem_obj = json.load(open("../../Objects/Base_Biochemistry.json"))
    biochem_obj['compounds'] = parse_compounds('../../Biochemistry/compounds.json')
    biochem_obj['reactions'] = parse_reactions('../../Biochemistry/reactions.json')

    biochem_obj['compound_aliases'] = defaultdict(lambda: defaultdict(list))
    for file in ("Compounds_Aliases.tsv","Names_Compounds_Aliases.tsv"):
        aliases = parse_aliases('../../Biochemistry/Aliases/'+file)
        for msid in aliases:
            for source in aliases[msid]:
                biochem_obj['compound_aliases'][msid][source].append(aliases[msid][source])
                
    biochem_obj['reaction_aliases'] = parse_aliases('../../Biochemistry/Aliases/Reactions_Aliases.tsv')
    biochem_obj['description'] = \
        'Biochemistry object generated on {} from ModelSEEDDatabase commit {}'\
        .format(datetime.today(), local_head)
    out_fp = '../../Objects/{}.json'.format(biochem_obj['id'])
    json.dump(biochem_obj, open(out_fp, 'w'), indent=4, sort_keys=True)
