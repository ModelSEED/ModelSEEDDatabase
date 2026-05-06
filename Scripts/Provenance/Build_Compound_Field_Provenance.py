#!/usr/bin/env python
"""
Build per-shard compound_NN.provenance.tsv sidecar files capturing per-field
source attribution for every ModelSEED compound.

Output columns (per compound):
  id                Compound ID
  formula_src       Source(s) of the picked structure that supplied formula
  charge_src        Source(s) of the picked structure that supplied charge
  inchikey_src      Source(s) of the picked InChIKey
  smiles_src        Source(s) of the picked SMILES
  pka_src           Source whose pKa entry was applied
  pkb_src           Source whose pKb entry was applied
  deltag_method     Which method's energy was used as deltag: GC, EQ, or empty
  deltag_used_src   Provenance pointer for deltag_method
  deltag_gc_src     Provenance pointer for thermodynamics.GroupContribution
  deltag_eq_src     Provenance pointer for thermodynamics.eQuilibrator

Provenance entries are formatted as:
  <DB>:<external_id>            when the picked structure is the same in both
                                Charged and Original protonation stages
  <DB>:<external_id>@<stage>    when only one stage produced the structure
  Override:<reason>             when a manual curation override was applied
  eQuilibrator:<MNXM_id>        for eQuilibrator-derived energies

The script does not recompute anything: it joins existing files and replays
the attribution logic of the four production update scripts:
  Scripts/Structures/Update_Compound_Structures_Formulas_Charge.py
  Scripts/Structures/Update_Compound_pKas.py
  Scripts/Thermodynamics/Update_Compound_GroupContribution_Energies.py
  Scripts/Thermodynamics/Update_Compound_eQuilibrator_Energies.py
"""
import os
import sys
import glob

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(SCRIPT_DIR, '..', '..', 'Libs', 'Python'))
from BiochemPy import Compounds  # noqa: E402

BIOCHEM_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', 'Biochemistry'))
STRUCT_ROOT  = os.path.join(BIOCHEM_ROOT, 'Structures')
THERMO_ROOT  = os.path.join(BIOCHEM_ROOT, 'Thermodynamics')
ACPS_FILE    = os.path.normpath(os.path.join(BIOCHEM_ROOT, 'Curation', 'overrides', 'acps_formula_charge.tsv'))
IGNORES_DIR  = os.path.normpath(os.path.join(BIOCHEM_ROOT, 'Curation', 'ignores'))

STRUCTURE_SOURCES   = ['KEGG', 'MetaCyc', 'ChEBI', 'Rhea']
PRIORITY_ORDER      = ['MetaCyc', 'KEGG', 'ChEBI', 'Rhea']
PKA_REPLAY_SOURCES  = ['KEGG', 'MetaCyc']  # matches Update_Compound_pKas.py hardcode
GC_REPLAY_SOURCES   = ['KEGG', 'MetaCyc']
EXCLUDED_CPDS       = {'cpd11632'}  # 'Light', hardcoded skip in production

HEADERS = [
    'id',
    'formula_src', 'charge_src', 'inchikey_src', 'smiles_src',
    'pka_src', 'pkb_src',
    'deltag_method', 'deltag_used_src', 'deltag_gc_src', 'deltag_eq_src',
]


def fmt(source, ext_id, stage=None):
    if stage is None:
        return f"{source}:{ext_id}"
    return f"{source}:{ext_id}@{stage}"


def load_all_evidence(path):
    """All_ModelSEED_Structures.txt → evidence[cpd][type][structure] = [(source, stage, ext_id), ...]"""
    evidence = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            cpd, stype, stage, ext_id, source, _formula, _charge, structure = parts[:8]
            evidence.setdefault(cpd, {}).setdefault(stype, {}).setdefault(structure, []).append((source, stage, ext_id))
    return evidence


def evidence_for(all_evidence, cpd, stype, structure):
    """Format evidence for a picked structure as ';'-joined DB:ext@stage entries."""
    if structure in (None, '', 'null'):
        return ''
    rows = all_evidence.get(cpd, {}).get(stype, {}).get(structure, [])
    if not rows:
        return ''
    grouped = {}
    for source, stage, ext_id in rows:
        grouped.setdefault((source, ext_id), set()).add(stage)
    priority_idx = {s: i for i, s in enumerate(PRIORITY_ORDER)}
    items = sorted(grouped.items(), key=lambda x: (priority_idx.get(x[0][0], 99), x[0][0], x[0][1]))
    out = []
    for (source, ext_id), stages in items:
        if stages == {'Charged', 'Original'}:
            out.append(fmt(source, ext_id))
        else:
            stage = next(iter(stages))
            out.append(fmt(source, ext_id, stage))
    return ';'.join(out)


def load_pka_data(helper, sources):
    """{(DB, ext_id): {'pKa': value, 'pKb': value}}.

    Reads from the post-A1 layout via Compounds.loadPerSourcePkas, which
    pulls every <db>/pkas/*.tsv. Returns the same dict shape as before
    so the downstream replay logic is unchanged.
    """
    data = helper.loadPerSourcePkas(sources)
    seen = sorted({db for (db, _ext) in data.keys()})
    return data, seen


def load_gc_data(thermo_root, sources):
    """{ext_id: (source, stage, dg_str)} matching the production load order:
       KEGG_Charged → KEGG_Original → MetaCyc_Charged → MetaCyc_Original.
       Charged wins; Original fills in only if Charged was the 1e+07 default."""
    data = {}
    for source in sources:
        for stage in ['Charged', 'Original']:
            path = os.path.join(thermo_root, 'ModelSEED', f'{source}_{stage}_MolAnalysis.tbl')
            if not os.path.isfile(path):
                continue
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) < 9:
                        continue
                    ext_id = parts[0]
                    dg     = parts[7]
                    if ext_id not in data:
                        data[ext_id] = (source, stage, dg)
                    else:
                        prev = data[ext_id][2]
                        try:
                            prev_is_default = abs(float(prev) - 10000000.0) < 1.0
                        except ValueError:
                            prev_is_default = False
                        if prev_is_default and dg != '1e+07':
                            data[ext_id] = (source, stage, dg)
    return data


def load_mnx_map(struct_root):
    """{deprotonated_inchikey: [mnx_id, ...]}"""
    mp = {}
    path = os.path.join(struct_root, 'MetaNetX', 'Structures_in_ModelSEED_and_eQuilibrator.txt')
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            mnx, ikey = parts[0], parts[1]
            depro = '-'.join(ikey.split('-')[:2])
            mp.setdefault(depro, []).append(mnx)
    return mp


def load_eq_energies(thermo_root):
    """{mnx_id: dg_float}"""
    energies = {}
    path = os.path.join(thermo_root, 'eQuilibrator', 'MetaNetX_Compound_Energies.tbl')
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            mnx, dg, _dge = parts[0], parts[1], parts[2]
            if 'energy' in dg or dg == 'nan':
                continue
            try:
                energies[mnx] = float(dg)
            except ValueError:
                continue
    return energies


def load_acp_overrides(path):
    """ACP override file → set of cpd ids and the fields each overrides."""
    overrides = {}
    if not os.path.isfile(path):
        return overrides
    with open(path) as fh:
        header = None
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if header is None:
                header = parts
                continue
            cpd = parts[0]
            cols = {}
            for i in range(1, len(parts)):
                if i < len(header) and parts[i] not in ('null', '10000000', ''):
                    cols[header[i]] = parts[i]
            if cols:
                overrides[cpd] = cols
    return overrides


def load_ignored_cpds(ignores_dir, helper):
    """Cpds whose missing-structure rows get the R/0 override.
    Replays Update_Compound_Structures_Formulas_Charge.py: only
    Ignored_Structures_Publication2020.txt with Accepted=='None', mapped via
    MetaCyc aliases."""
    ignored = set()
    path = os.path.join(ignores_dir, 'Ignored_Structures_Publication2020.txt')
    if not os.path.isfile(path):
        return ignored
    smile_dict = helper.loadStructures(['SMILE', 'InChIKey'], ['KEGG', 'MetaCyc']).get('SMILE', {})
    src_aliases = helper.loadSourceAliases().get('MetaCyc', {})
    with open(path) as fh:
        first = True
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if first:
                first = False
                if parts[0] == 'ID':
                    continue
            if len(parts) < 2 or parts[1] != 'None':
                continue
            ext_id = parts[0]
            if ext_id in smile_dict:
                cpds = src_aliases.get(ext_id, [])
                if len(cpds) == 1:
                    ignored.add(cpds[0])
    return ignored


def notes_set(notes):
    if notes in (None, '', 'null'):
        return set()
    if isinstance(notes, list):
        return {n for n in notes if n}
    if isinstance(notes, str):
        return {n for n in notes.split('|') if n}
    return set()


def main():
    print('Loading compounds and aliases...')
    helper    = Compounds()
    compounds = helper.loadCompounds()
    aliases   = helper.loadMSAliases()

    print('Loading All_ModelSEED_Structures.txt...')
    all_evidence = load_all_evidence(os.path.join(STRUCT_ROOT, 'All_ModelSEED_Structures.txt'))

    print('Loading Unique_ModelSEED_Structures.txt picks...')
    structures = helper.loadStructures(['SMILE', 'InChI', 'InChIKey'], ['ModelSEED'])

    print('Loading curation overrides...')
    acp_overrides = load_acp_overrides(ACPS_FILE)
    ignored_cpds  = load_ignored_cpds(IGNORES_DIR, helper)

    print('Loading pKa data...')
    pka_data, pka_sources_seen = load_pka_data(helper, STRUCTURE_SOURCES)
    print(f'  pKa files present: {pka_sources_seen}')
    print(f'  pKa replay sources (matching production): {PKA_REPLAY_SOURCES}')

    print('Loading thermodynamics (GC, eQ)...')
    gc_data     = load_gc_data(THERMO_ROOT, GC_REPLAY_SOURCES)
    mnx_by_ikey = load_mnx_map(STRUCT_ROOT)
    eq_energies = load_eq_energies(THERMO_ROOT)

    print(f'Building provenance for {len(compounds)} compounds...')

    provenance = {}
    for cpd in sorted(compounds.keys()):
        cpd_obj = compounds[cpd]
        prov = {k: '' for k in HEADERS[1:]}

        if cpd in EXCLUDED_CPDS:
            provenance[cpd] = prov
            continue

        # ---- formula / charge: pick InChI first key, fallback to SMILE first key ----
        struct_picked_type = None
        struct_picked      = None
        if cpd in structures:
            if 'InChI' in structures[cpd]:
                struct_picked_type = 'InChI'
                struct_picked      = list(structures[cpd]['InChI'].keys())[0]
            elif 'SMILE' in structures[cpd]:
                struct_picked_type = 'SMILE'
                struct_picked      = list(structures[cpd]['SMILE'].keys())[0]

        if struct_picked is not None:
            ev = evidence_for(all_evidence, cpd, struct_picked_type, struct_picked)
            prov['formula_src'] = ev
            prov['charge_src']  = ev
        elif cpd in ignored_cpds:
            prov['formula_src'] = 'Override:ignored_R_group'
            prov['charge_src']  = 'Override:ignored_R_group'

        # ACP override is layered on top of the structure-derived value
        if cpd in acp_overrides:
            ovr = acp_overrides[cpd]
            if 'formula' in ovr:
                prov['formula_src'] = (prov['formula_src'] + ';' if prov['formula_src'] else '') + 'Override:ACPs'
            if 'charge' in ovr:
                prov['charge_src']  = (prov['charge_src']  + ';' if prov['charge_src']  else '') + 'Override:ACPs'

        # ---- inchikey: sorted-first across stages ----
        if cpd in structures and 'InChIKey' in structures[cpd]:
            ikey = sorted(structures[cpd]['InChIKey'].keys())[0]
            prov['inchikey_src'] = evidence_for(all_evidence, cpd, 'InChIKey', ikey)

        # ---- smiles: sorted-first across stages ----
        if cpd in structures and 'SMILE' in structures[cpd]:
            smi = sorted(structures[cpd]['SMILE'].keys())[0]
            prov['smiles_src'] = evidence_for(all_evidence, cpd, 'SMILE', smi)

        # ---- pKa / pKb: KEGG-then-MetaCyc, first alias with a hit ----
        for db in PKA_REPLAY_SOURCES:
            if db not in aliases.get(cpd, {}):
                continue
            for ext_id in aliases[cpd][db]:
                if (db, ext_id) in pka_data:
                    entry = pka_data[(db, ext_id)]
                    if 'pKa' in entry and not prov['pka_src']:
                        prov['pka_src'] = fmt(db, ext_id)
                    if 'pKb' in entry and not prov['pkb_src']:
                        prov['pkb_src'] = fmt(db, ext_id)
                    break
            if prov['pka_src'] or prov['pkb_src']:
                break

        # ---- GC ΔG: replay alias-filter against picked structure, lowest dg wins ----
        gc_struct_type = None
        if cpd in structures:
            if 'InChIKey' in structures[cpd]:
                gc_struct_type = 'InChIKey'
            elif 'SMILE' in structures[cpd]:
                gc_struct_type = 'SMILE'
        if gc_struct_type:
            struct = list(structures[cpd][gc_struct_type].keys())[0]
            valid_aliases = sorted({ext for _src, _stg, ext in all_evidence.get(cpd, {}).get(gc_struct_type, {}).get(struct, [])})
            lowest_dg  = None
            lowest_src = None
            for ext_id in valid_aliases:
                if ext_id not in gc_data:
                    continue
                src, stg, dg = gc_data[ext_id]
                try:
                    dg_f = float(dg)
                except ValueError:
                    continue
                if lowest_dg is None or dg_f < lowest_dg:
                    lowest_dg  = dg_f
                    lowest_src = (src, stg, ext_id)
            if lowest_src:
                src, stg, ext_id = lowest_src
                prov['deltag_gc_src'] = fmt(src, ext_id, stg)

        # ---- eQ ΔG: cpd → first sorted InChIKey/SMILE → first 2 layers → MNX → lowest energy ----
        eq_struct_type = None
        if cpd in structures:
            if 'InChIKey' in structures[cpd]:
                eq_struct_type = 'InChIKey'
            elif 'SMILE' in structures[cpd]:
                eq_struct_type = 'SMILE'
        if eq_struct_type:
            struct = list(structures[cpd][eq_struct_type].keys())[0]
            depro  = '-'.join(struct.split('-')[:2])
            mnxs   = mnx_by_ikey.get(depro, [])
            lowest_dg  = None
            lowest_mnx = None
            for mnx in mnxs:
                if mnx not in eq_energies:
                    continue
                dg = eq_energies[mnx]
                if lowest_dg is None or dg < lowest_dg:
                    lowest_dg  = dg
                    lowest_mnx = mnx
            if lowest_mnx:
                prov['deltag_eq_src'] = fmt('eQuilibrator', lowest_mnx)

        # ---- ΔG method: read existing notes (GC|EQ|EQU encoding) ----
        n = notes_set(cpd_obj.get('notes'))
        if 'EQU' in n:
            prov['deltag_method']   = 'EQ'
            prov['deltag_used_src'] = prov['deltag_eq_src']
        elif 'GC' in n:
            prov['deltag_method']   = 'GC'
            prov['deltag_used_src'] = prov['deltag_gc_src']
        # else: leave both empty

        provenance[cpd] = prov

    # ---- Write per-shard files matching saveCompounds bucketing ----
    print('Writing provenance shards...')
    cpds_root = os.path.join(BIOCHEM_ROOT, 'compound_')
    cpds_count_thousands = 0
    fh = open(f'{cpds_root}{cpds_count_thousands:02}.provenance.tsv', 'w')
    fh.write('\t'.join(HEADERS) + '\n')
    prev_rounded = 0

    for cpd in sorted(provenance.keys()):
        cpd_count    = int(cpd[3:])
        cur_rounded  = cpd_count - cpd_count % 1000
        if cur_rounded > prev_rounded:
            prev_rounded = cur_rounded
            cpds_count_thousands += 1
            fh.close()
            fh = open(f'{cpds_root}{cpds_count_thousands:02}.provenance.tsv', 'w')
            fh.write('\t'.join(HEADERS) + '\n')
        prov = provenance[cpd]
        row  = [cpd] + [prov[k] for k in HEADERS[1:]]
        fh.write('\t'.join(row) + '\n')
    fh.close()
    print('Done.')


if __name__ == '__main__':
    main()
