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
        'abbreviation': ('abbreviation', str, ''), #checked v1.0
        'deltag': ('deltaG', float, 10000000), #checked v1.0
        'deltagerr': ('deltaGErr', float, 10000000), #checked v1.0
        'formula': ('formula', str, ''), #checked v1.0
        'id': ('id', str, ''), #checked v1.0
        'name': ('name', str, ''), #checked v1.0
        'mass': ('mass', float, 0), #checked v1.0
        'charge': ('defaultCharge', int, 10000000), #checked v1.0
        'is_cofactor': ('isCofactor', int, 0), #checked v1.0
    }
    source = json.load(open(fp))
    compounds = []
    for in_comp in source:
#        if in_comp['is_obsolete']:
#            continue
        out_comp = {"unchargedFormula": "", 'cues': {}}
        for in_id in trans_dict.keys():
            try:
                out_comp[trans_dict[in_id][0]] = trans_dict[in_id][1](in_comp[in_id])
            except ValueError:
                print(str(in_comp[in_id])+" is right type, but wrong value, for "+in_id+" in "+in_comp['id'])
                out_comp[trans_dict[in_id][0]] = trans_dict[in_id][1](trans_dict[in_id][2])
            except TypeError:
                print(str(in_comp[in_id])+" is wrong type for "+in_id+" in "+in_comp['id'])
                out_comp[trans_dict[in_id][0]] = trans_dict[in_id][1](trans_dict[in_id][2])

        #NB: pKa/pKb only stored for first molecular fragment
        out_comp['pkas'] = dict([(int(x.split(':')[1]), [float(x.split(':')[2])])
                                 for x in in_comp['pka'].split(';') if x and x[0] == '1'])
        out_comp['pkbs'] = dict([(int(x.split(':')[1]), [float(x.split(':')[2])])
                                 for x in in_comp['pkb'].split(';') if x and x[0] == '1'])
        compounds.append(out_comp)
    return compounds

#{'reagents': [{'isCofactor': 0, 'compound_ref': '~/compounds/id/cpd12036', 'compartment_ref': '~/compartments/id/c', 'coefficient': 1}, {'isCofactor': 0, 'compound_ref': '~/compounds/id/cpd11463', 'compartment_ref': '~/compartments/id/c', 'coefficient': -1}], 'deltaG': 10000000, 'thermoReversibility': '=', 'status': 'MI:C12/N7/O1/R-1|CI:3|HI:-21', 'name': 'rxn05785', 'cues': {'~/cues/id/WOCOW': 1, '~/cues/id/RWCHdblW': 2, '~/cues/id/WdblNH2': 1, '~/cues/id/WWCHW': -1, '~/cues/id/WdblCWW': 1, '~/cues/id/WNH3': 1, '~/cues/id/RWOW': 1, '~/cues/id/HeteroAromatic': 2, '~/cues/id/RWdblNW': 3, '~/cues/id/RWCHWW': 4, '~/cues/id/amide': -1, '~/cues/id/WCOOn': -1, '~/cues/id/RWWNW': 1, '~/cues/id/PrimOH': 1, '~/cues/id/RWCdblWW': 1, '~/cues/id/TWWCdblW': 2, '~/cues/id/WCH2W': 4, '~/cues/id/WketoneW': -1, '~/cues/id/RW': -1, '~/cues/id/WNH2': 1}, 'direction': '=', 'abbreviation': '', 'id': 'rxn05785', 'deltaGErr': 10000000}

def parse_reactions(fp):
    # mapping of field in old to field in new and type
    trans_dict = {
        'abbreviation': ('abbreviation', str, ''),
        'deltag': ('deltaG', float, 10000000),
        'deltagerr': ('deltaGErr', float, 10000000),
        'direction': ('direction', str, '='),
        'id': ('id', str, ''),
        'name': ('name', str, ''),
        'reversibility': ('thermoReversibility', str, '='),
        'status': ('status', str, ''),
    }
    source = json.load(open(fp))
    reactions = []
    for in_rxn in source:
#        if in_rxn['is_obsolete']:
#            continue
        out_rxn = {"defaultProtons": 0, 'cues': {}}
        for in_id in trans_dict.keys():
            try:
                out_rxn[trans_dict[in_id][0]] = trans_dict[in_id][1](in_rxn[in_id])
            except ValueError:
                print(str(in_rxn[in_id])+" is right type, but wrong value, for "+in_id+" in "+in_rxn['id'])
                out_rxn[trans_dict[in_id][0]] = trans_dict[in_id][1](trans_dict[in_id][2])
            except TypeError:
                print(str(in_rxn[in_id])+" is wrong type for "+in_id+" in "+in_rxn['id'])
                out_rxn[trans_dict[in_id][0]] = trans_dict[in_id][1](trans_dict[in_id][2])

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
    if r.fieldnames != ['ModelSEED ID', 'External ID', 'Source']:
        raise ValueError("Invalid aliases table format")
    aliases = defaultdict(lambda: defaultdict(list))
    for line in r:
        aliases[line['ModelSEED ID']][line['Source']].append(line['External ID'])
    return aliases


if __name__ == "__main__":
    biochem_obj = json.load(open("../../Objects/Base_Biochemistry.json"))
    biochem_obj['compounds'] = parse_compounds('../../Biochemistry/compounds.json')
    biochem_obj['reactions'] = parse_reactions('../../Biochemistry/reactions.json')

    biochem_obj['compound_aliases'] = defaultdict(lambda: defaultdict(list))
    for file in ("Unique_ModelSEED_Compound_Aliases.txt","Unique_ModelSEED_Compound_Names.txt"):
        aliases = parse_aliases('../../Biochemistry/Aliases/'+file)
        for msid in aliases:
            for source in aliases[msid]:
                biochem_obj['compound_aliases'][msid][source]=aliases[msid][source]
                
    biochem_obj['reaction_aliases'] = parse_aliases('../../Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt')
    biochem_obj['description'] = \
        'Biochemistry object generated on {} from ModelSEEDDatabase commit {}'\
        .format(datetime.today(), local_head)
    out_fp = '../../Objects/{}.json'.format(biochem_obj['id']+"_Biochem")
    json.dump(biochem_obj, open(out_fp, 'w'), indent=4, sort_keys=True)
