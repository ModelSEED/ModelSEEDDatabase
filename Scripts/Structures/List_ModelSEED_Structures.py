#!/usr/bin/env python
import os
import sys
import csv
import json
import glob

#################################################################
## Load Compound Objects into memory
#################################################################
sys.path.append('../../Libs/Python')
from BiochemPy import Compounds


#################################################################
## Curated structure picks (manual overrides)
#################################################################

def load_curated_picks():
    """Read all *.tsv files under Biochemistry/Curation/overrides/structure_picks/
    and return {cpd_id: {format, structure, source_db, source_id, curator,
                          date, rationale}}.

    Each TSV file's name is the curator (e.g. samseaver.tsv). Last-write-
    wins if a compound appears in multiple curator files; warns on overlap.
    """
    picks = {}
    overrides_dir = os.path.join(
        os.path.dirname(__file__), '..', '..',
        'Biochemistry', 'Curation', 'overrides', 'structure_picks')
    overrides_dir = os.path.normpath(overrides_dir)
    if not os.path.isdir(overrides_dir):
        return picks
    for path in sorted(glob.glob(os.path.join(overrides_dir, '*.tsv'))):
        curator = os.path.splitext(os.path.basename(path))[0]
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                cpd_id = row.get('cpd_id', '').strip()
                if not cpd_id:
                    continue
                if cpd_id in picks and picks[cpd_id]['curator'] != curator:
                    print(f"WARN: cpd {cpd_id} curated by both "
                          f"{picks[cpd_id]['curator']} and {curator}; "
                          f"using {curator}", file=sys.stderr)
                picks[cpd_id] = {
                    'format': row.get('format', '').strip(),
                    'structure': row.get('structure', '').strip(),
                    'source_db': row.get('source_db', '').strip(),
                    'source_id': row.get('source_id', '').strip(),
                    'curator': curator,
                    'date': row.get('date', '').strip(),
                    'rationale': row.get('rationale', '').strip(),
                }
    return picks


def derive_structures_from_override(override):
    """Given a curator override row, parse its structure via RDKit and
    return {'SMILE': str, 'InChI': str, 'InChIKey': str, 'formula': str,
             'charge': str} -- all derived consistently from the override.
    Returns None if RDKit can't parse the structure.
    """
    import re
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdMolDescriptors
        from rdkit.Chem.inchi import MolFromInchi, MolToInchi, InchiToInchiKey
        RDLogger.DisableLog('rdApp.*')
    except ImportError:
        print("ERROR: RDKit not available; cannot derive structures "
              "from curator overrides. Skipping override consult.",
              file=sys.stderr)
        return None

    fmt = override['format']
    struct = override['structure']
    if fmt == 'InChI':
        mol = MolFromInchi(struct)
        if mol is None:
            return None
        inchi = struct
    elif fmt == 'SMILE':
        mol = Chem.MolFromSmiles(struct)
        if mol is None:
            return None
        # Standard InChI can't represent wildcard atoms (*). For
        # R-containing picks we still return SMILES + formula + charge;
        # InChI and InChIKey are left blank and downstream skips them.
        try:
            inchi = MolToInchi(mol) or ''
        except Exception:
            inchi = ''
    else:
        return None

    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    inchikey = InchiToInchiKey(inchi) if inchi else ''
    # rdMolDescriptors.CalcMolFormula: strip trailing charge marker
    # (stored in a separate column), rewrite '*' wildcards as 'R'
    # (ModelSEED convention), and re-sort into Hill+alpha order --
    # RDKit places '*' immediately after H, but the rest of the
    # codebase places R alphabetically among the heavy atoms
    # (e.g. C22H31N7O17P3RS, not C22H31RN7O17P3S).
    formula_raw = re.sub(r'[+\-]\d*$', '',
                         rdMolDescriptors.CalcMolFormula(mol))
    counts = {}
    for m in re.finditer(r'([A-Z][a-z]?|\*)(\d*)', formula_raw):
        el = 'R' if m.group(1) == '*' else m.group(1)
        counts[el] = counts.get(el, 0) + int(m.group(2) or 1)
    def _key(el):
        if el == 'C': return (0, '')
        if el == 'H': return (1, '')
        return (2, el)
    formula = ''.join(
        f"{el}{counts[el] if counts[el] > 1 else ''}"
        for el in sorted(counts, key=_key))
    charge = str(Chem.GetFormalCharge(mol))
    return {
        'SMILE': smiles,
        'InChI': inchi,
        'InChIKey': inchikey,
        'formula': formula,
        'charge': charge,
    }


CURATED_PICKS = load_curated_picks()


#################################################################
## ACP formula/charge overrides
##
## Biochemistry/Curation/overrides/acps_formula_charge.tsv contains
## hand-curated formula/charge values for ACP (acyl-carrier-protein)
## compounds. The formulas include the pantetheine + 4'-phosphate side
## chain plus the specific acyl group, with a wildcard 'R' representing
## the protein backbone. Historically only Update_Compound_Structures_
## Formulas_Charge.py consulted this file (to override the per-compound
## JSON records); the picker did not, so ACP compounds hit the
## formula-conflict branch and were dropped from Unique_ModelSEED_
## Structures.txt entirely (including cpd11493 ACP itself, used in
## 911 reactions).
##
## The picker now reads this file too and uses it as the source of
## truth for formula/charge on any compound listed. The picked
## structure comes from a source database (typically KEGG or MetaCyc)
## because the override file records formula/charge only, not the
## actual SMILES/InChI/InChIKey.
#################################################################

def load_acps_overrides():
    """Read Biochemistry/Curation/overrides/acps_formula_charge.tsv
    and return {cpd_id: (formula, charge)}."""
    path = os.path.normpath(os.path.join(
        os.path.dirname(__file__), '..', '..',
        'Biochemistry', 'Curation', 'overrides', 'acps_formula_charge.tsv'))
    d = {}
    if not os.path.isfile(path):
        return d
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            cpd = row.get('ID', '').strip()
            if not cpd:
                continue
            d[cpd] = (row.get('formula', '').strip(),
                      row.get('charge', '').strip())
    return d


ACP_OVERRIDES = load_acps_overrides()


#################################################################
## SMILES canonicalization (idempotent with Recanonicalize_SMILES.py)
##
## Every SMILES row this script writes -- to All_ModelSEED_Structures.txt
## and to Unique_ModelSEED_Structures.txt -- is passed through the same
## RDKit-canonical function that Recanonicalize_SMILES.py uses. This
## makes the two scripts' outputs idempotent: running the picker followed
## by Recanonicalize produces zero further changes.
##
## Import lazily-guarded so the picker still runs (with warnings) if
## RDKit isn't installed; the SMILES would then be written as-is.
#################################################################

def _make_smiles_canonicalizer():
    """Return a function s -> canonical(s) that never raises;
    on parse failure it returns the input unchanged."""
    try:
        # Reuse the algorithm from Recanonicalize_SMILES.py in the same
        # directory to guarantee identical output.
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if script_dir not in sys.path:
            sys.path.insert(0, script_dir)
        from Recanonicalize_SMILES import canonical_smiles as _canon
        def _fn(s):
            if not s:
                return s
            new, status = _canon(s)
            return new if status not in ('parse_failure', 'null_or_empty') else s
        return _fn
    except Exception as e:
        print(f"WARN: could not load SMILES canonicalizer ({e}); "
              f"picker will write SMILES as-is from sources.", file=sys.stderr)
        return lambda s: s


CANONICALIZE_SMILES = _make_smiles_canonicalizer()

#Load Compounds
CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()

#################################################################
## Load Formula Strings from file
#################################################################

Structures_Root=os.path.dirname(__file__)+"/../../Biochemistry/Structures/"
Formulas_Dict = CompoundsHelper.loadPerSourceFormulasCharges(['InChI','SMILE'], ['KEGG','MetaCyc','ChEBI','Rhea'])

#################################################################
## Load Curated Picks for Structures
#################################################################

#Load Curated Structures
Ignored_Structures=dict()
for file in glob.glob(os.path.dirname(__file__)+"/../../Biochemistry/Curation/ignores/*.txt"):
    with open(file) as ignore_file:
        for line in ignore_file.readlines():
            array=line.split('\t')
            Ignored_Structures[array[0]]=1

#################################################################
## Load Structures and Aliases
#################################################################

#Load Structures and Aliases
Structures_Dict = CompoundsHelper.loadStructures(["SMILE","InChIKey","InChI"],["KEGG","MetaCyc","ChEBI","Rhea"])
MS_Aliases_Dict =  CompoundsHelper.loadMSAliases(["KEGG","MetaCyc","ChEBI","Rhea"])

#################################################################
## Open filehandles for writing
#################################################################

master_structs_file = open(Structures_Root+"All_ModelSEED_Structures.txt",'w')
unique_structs_file = open(Structures_Root+"Unique_ModelSEED_Structures.txt",'w')
unique_structs_file.write("ID\tType\tAliases\tFormula\tCharge\tStructure\n")
structure_conflicts_file = open(os.path.dirname(__file__)+"/../../Biochemistry/Structures/_reports/Structure_Conflicts.txt",'w')
formula_conflicts_file = open(os.path.dirname(__file__)+"/../../Biochemistry/Structures/_reports/Formula_Conflicts.txt",'w')
pick_reasons_file = open(Structures_Root+"Pick_Reasons.txt",'w')
pick_reasons_file.write("ID\tType\tStage\tReason\tChosen_Structure\tChosen_Aliases\n")

#################################################################
## Iterate through ModelSEED identifiers
#################################################################

for msid in sorted(MS_Aliases_Dict.keys()):

    #################################################################
    ## For the ModelSEED compound build dict of all structures
    #################################################################

    Structs = dict()
    Formulas=dict()
    for source in 'KEGG','MetaCyc','ChEBI','Rhea':
        if(source not in MS_Aliases_Dict[msid].keys()):
            continue

        #################################################################
        ## Iterate through types, sources, ids
        #################################################################

        for struct_type in sorted(Structures_Dict.keys()):
            for external_id in sorted(MS_Aliases_Dict[msid][source]):
                if(external_id not in Structures_Dict[struct_type]):
                    continue

                for struct_stage in sorted(Structures_Dict[struct_type][external_id].keys()):
                    if(struct_type not in Structs):
                        Structs[struct_type]=dict()
                        Formulas[struct_type]=dict()

                    if(struct_stage not in Structs[struct_type]):
                        Structs[struct_type][struct_stage]=dict()
                        Formulas[struct_type][struct_stage]=dict()

                    for structure in sorted(Structures_Dict[struct_type][external_id][struct_stage].keys()):

                        formula_charge_dict={'formula':"null",'charge':"null"}

                        if(struct_type in Formulas_Dict[source] and external_id in Formulas_Dict[source][struct_type][struct_stage]):
                            formula_charge_dict = Formulas_Dict[source][struct_type][struct_stage][external_id]

                        #################################################################
                        ## Write all structures to master 'All_ModelSEED_Strutures.txt'
                        ## SMILES rows are canonicalized so this file's output stays
                        ## idempotent with Recanonicalize_SMILES.py.
                        #################################################################

                        structure_to_write = (CANONICALIZE_SMILES(structure)
                                              if struct_type == "SMILE" else structure)
                        master_structs_file.write("\t".join([msid,struct_type,struct_stage,external_id,source,\
                                                                 formula_charge_dict['formula'],\
                                                                 formula_charge_dict['charge'],\
                                                                 structure_to_write])+"\n")
                        
                        #################################################################
                        ## Skip if curated structure designated to be ignored
                        #################################################################

                        if(external_id in Ignored_Structures):
                            continue

                        #################################################################
                        ## Populate structures dictionary
                        #################################################################

                        if(structure not in Structs[struct_type][struct_stage]):
                            Structs[struct_type][struct_stage][structure]=dict()
                        Structs[struct_type][struct_stage][structure][external_id]=source

                        formula_charge_json = json.dumps(formula_charge_dict)
                        if(formula_charge_json not in Formulas[struct_type][struct_stage]):
                            Formulas[struct_type][struct_stage][formula_charge_json]=dict()
                        Formulas[struct_type][struct_stage][formula_charge_json][external_id]=source
                        
    #################################################################
    ## Skip if no structures were collected
    #################################################################

    if(len(Structs.keys())==0):
        continue

    #################################################################
    ## Curator override consult: if this compound has a manual pick in
    ## Biochemistry/Curation/overrides/structure_picks/<curator>.tsv,
    ## use it directly and skip the cascade tiebreaker entirely.
    ## All three formats (SMILE, InChIKey, InChI) are derived from the
    ## override's structure via RDKit so they stay internally consistent.
    #################################################################

    if(msid in CURATED_PICKS):
        override = CURATED_PICKS[msid]
        derived = derive_structures_from_override(override)
        if(derived is None):
            print(f"WARN: cpd {msid} curator override ({override['curator']}) "
                  f"could not be parsed by RDKit; falling back to cascade",
                  file=sys.stderr)
        else:
            # Aggregate all aliases for this compound across all sources
            override_aliases = set()
            for src in MS_Aliases_Dict[msid]:
                for alias in MS_Aliases_Dict[msid][src]:
                    override_aliases.add(alias)
            aliases_str = ";".join(sorted(override_aliases))
            # Write the structure rows to Unique using derived values.
            # Skip any format the pick doesn't yield (e.g. InChI/InChIKey
            # for R-containing SMILES picks that RDKit can't represent).
            for stype in ("SMILE", "InChIKey", "InChI"):
                if not derived.get(stype):
                    continue
                unique_structs_file.write("\t".join((
                    msid, stype, aliases_str,
                    derived['formula'], derived['charge'],
                    derived[stype])) + "\n")
            # Record the pick rationale
            reason = f"manual_curation:{override['curator']}"
            pick_reasons_file.write("\t".join((
                msid, override['format'], "Charged", reason,
                override['structure'], aliases_str)) + "\n")
            # Also report the underlying conflict (transparency) if any
            # We pick the priority type/stage for the conflict report only.
            for try_type in ("InChI", "SMILE"):
                if(try_type in Structs and "Charged" in Structs[try_type]
                        and len(Structs[try_type]["Charged"]) > 1):
                    for structure in Structs[try_type]["Charged"]:
                        for ext_id in Structs[try_type]["Charged"][structure]:
                            structure_conflicts_file.write("\t".join((
                                msid, try_type, "Charged", structure, ext_id,
                                Structs[try_type]["Charged"][structure][ext_id])) + "\n")
                    break
            # Skip the cascade tiebreaker for this compound
            continue

    #################################################################
    ## ACP formula-override consult: for any compound with a hand-
    ## curated formula/charge in acps_formula_charge.tsv, apply the
    ## override and pick any source structure. This resolves compounds
    ## that would otherwise hit the formula_conflict_no_pick branch
    ## and be dropped from Unique.
    #################################################################

    if(msid in ACP_OVERRIDES):
        override_formula, override_charge = ACP_OVERRIDES[msid]
        # Prefer InChI/Charged, then InChI/Original, then SMILE/Charged, SMILE/Original
        pick_type, pick_stage = None, None
        for try_type in ('InChI', 'SMILE'):
            if try_type not in Structs:
                continue
            for try_stage in ('Charged', 'Original'):
                if try_stage in Structs[try_type]:
                    pick_type, pick_stage = try_type, try_stage
                    break
            if pick_type:
                break
        if pick_type is None:
            print(f"Warning: ACP override for {msid} but no source "
                  f"structure available; skipping.", file=sys.stderr)
        else:
            # Aggregate all aliases for this compound across sources
            override_aliases = set()
            for src in MS_Aliases_Dict[msid]:
                for alias in MS_Aliases_Dict[msid][src]:
                    override_aliases.add(alias)
            aliases_str = ";".join(sorted(override_aliases))
            # Write each format's Unique row using the picked source's structure
            # for that format, plus the override formula/charge.
            representative_structure = sorted(Structs[pick_type][pick_stage].keys())[0]
            for stype in ('SMILE', 'InChIKey', 'InChI'):
                if stype not in Structs or pick_stage not in Structs[stype]:
                    continue
                s = sorted(Structs[stype][pick_stage].keys())[0]
                s_out = CANONICALIZE_SMILES(s) if stype == 'SMILE' else s
                unique_structs_file.write("\t".join((
                    msid, stype, aliases_str,
                    override_formula, override_charge, s_out)) + "\n")
            pick_reasons_file.write("\t".join((
                msid, pick_type, pick_stage, "manual_formula_override:acps",
                representative_structure, aliases_str)) + "\n")
            # Report underlying conflict (transparency)
            for try_type in ("InChI", "SMILE"):
                if(try_type in Structs and pick_stage in Structs[try_type]
                        and len(Structs[try_type][pick_stage]) > 1):
                    for structure in Structs[try_type][pick_stage]:
                        for ext_id in Structs[try_type][pick_stage][structure]:
                            structure_conflicts_file.write("\t".join((
                                msid, try_type, pick_stage, structure, ext_id,
                                Structs[try_type][pick_stage][structure][ext_id])) + "\n")
                    break
            continue

    #################################################################
    ## Prioritized which type and stage for the structure for comparison
    ## Priority Order is:
    ## 1) Charged InChI
    ## 2) Original InChI
    ## 3) Charged SMILE
    ## 4) Original SMILE
    #################################################################

    struct_type=None
    struct_stage=None
    if("InChI" in Structs):
        struct_type="InChI"
        if("Charged" in Structs[struct_type]):
            struct_stage="Charged"
        elif("Original" in Structs[struct_type]):
            struct_stage="Original"
    elif("SMILE" in Structs):
        struct_type="SMILE"
        if("Charged" in Structs[struct_type]):
            struct_stage="Charged"
        elif("Original" in Structs[struct_type]):
            struct_stage="Original"

    if(struct_type is None or struct_stage is None):
        #At time of writing, this doesn't happen
        print("Warning: no structures used for "+msid)
        continue
                            
    #################################################################
    ## Now the prioritized type and stage has been established
    ## We look to see whether or not they have the same structure
    ## from different sources
    ##
    ## If there's only one structure string, or there's multiple
    ## structure strings that have the same formula, the structure
    ## Is considered a 'pass' NB: This will be re-written
    ## To accomodate for differing protonation states being depicted
    ## as charges
    #################################################################

    struct_pass=0
    struct_conflict=0
    formula_conflict=0
    pick_reason=None
    if(len(Structs[struct_type][struct_stage].keys())==1):
        struct_pass=1
        pick_reason="single_structure"
    elif(len(Formulas[struct_type][struct_stage].keys())==1):
        struct_conflict=1
        struct_pass=1
    else:
        struct_conflict=1
        formula_conflict=1
        pass
                            
    #################################################################
    ## Here on, we consider the structure(s) in general for a
    ## ModelSEED compound to be unique to that compound though
    ## they may vary if there's still a conflict.
    ## The code in this section describes how we
    ## might 'pick' a structure.
    #################################################################

    if(struct_pass):

        #################################################################
        ## The formula is considered to be identical across structures
        ## But we've recently found this isn't always the case with some
        ## protonated alcohol groups
        #################################################################

        #Only one formula/charge combination possible here
        formula_charge_dict=json.loads(list(Formulas[struct_type][struct_stage].keys())[0])
        
        #################################################################
        ## If there are no structural conflicts then all the structures
        ## are identical and we pick one
        #################################################################

        if(struct_conflict == 0):

            #################################################################
            ## Iterate through the types for writing to file
            ## Keeping order consistent but it doesn't matter which
            #################################################################

            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                #################################################################
                ## Even though there is one standardized InChI structure
                ## There can be more than one SMILE so here we choose one 
                #################################################################

                structure = sorted(Structs[structure_type][struct_stage].keys())[0]

                #################################################################
                ## But we collect all aliases for that type and stage
                ## This means that the aliases for different SMILE strings
                ## will be grouped together under a single SMILE string
                #################################################################
                
                aliases=dict()
                for struct in Structs[structure_type][struct_stage].keys():
                    for alias in Structs[structure_type][struct_stage][struct]:
                        aliases[alias]=1

                #################################################################
                ## We write them to file (SMILES canonicalized to keep Unique
                ## idempotent with Recanonicalize_SMILES.py)
                #################################################################

                structure_to_write = (CANONICALIZE_SMILES(structure)
                                      if structure_type == "SMILE" else structure)
                unique_structs_file.write("\t".join((msid,\
                                                     structure_type,\
                                                     ";".join(sorted(aliases)),\
                                                     formula_charge_dict['formula'],\
                                                     formula_charge_dict['charge'],\
                                                     structure_to_write))+"\n")

        else:

            #################################################################
            ## Now here we have similar structures that are identical in
            ## formula but differ in connectivity somehow. The code in
            ## this section is about 'selecting' the primary structure from
            ## a database and an alias that represents the compound
            #################################################################
            
            struct_conflicts = dict()
            sources_structures=dict()

            #################################################################
            ## First we pick InChIKey over SMILE as the primary structure, if
            ## we can. Then we collect that type of structure
            #################################################################

            for structure_type in "InChIKey","SMILE":
                if(structure_type not in Structs):
                    continue

                for structure in Structs[structure_type][struct_stage]:
                    for alias in Structs[structure_type][struct_stage][structure]:
                        source = Structs[structure_type][struct_stage][structure][alias]
                        if(structure not in struct_conflicts):
                            struct_conflicts[structure]=dict()
                        if(source not in struct_conflicts):
                            struct_conflicts[structure][source]=dict()
                        struct_conflicts[structure][source][alias]=1
                        if(source not in sources_structures):
                            sources_structures[source]=dict()
                        sources_structures[source][structure]=1


                #################################################################
                ## Break if collected InChIKey structures
                #################################################################
                        
                if(len(struct_conflicts)>0):
                    break

            #################################################################
            ## Here we go through the structures and determine if any come
            ## from more than database to prioritize, and then if not, we
            ## break them down to see how their connectivity matches and then
            ## we prioritize certain databases.
            #################################################################

            chosen_structure = None

            #################################################################
            ## Find structures that are identical in more than one database
            #################################################################

            chosen_structures = dict()
            for structure in struct_conflicts:
                if(len(struct_conflicts[structure])>1):
                    chosen_structures[structure]=1

            if(len(chosen_structures)>0):
                if(len(chosen_structures)==1):

                        #################################################################
                        ## If there is one structure that is identical in more than one
                        ## database we pick that structure to be the primary one
                        #################################################################

                        chosen_structure = list(chosen_structures.keys())[0]
                        pick_reason="multi_source_agreement"

                else:
                    for structure in chosen_structures:

                        #################################################################
                        ## If more than one identical structures are found in more than
                        ## one database then we have to arbitrarily pick one
                        ## There's not many at all
                        #################################################################

                        #Avoid lack of stereochemistry, at time of writing, never happens
                        #For SMILE string
                        if('UHFFFAOYSA' not in structure):
                            chosen_structure = structure
                            pick_reason="multi_source_stereochemistry"
                            break

                    if(chosen_structure is None):

                        print(msid,chosen_structures)

                        chosen_structure = list(chosen_structures.keys())[0]
                        pick_reason="multi_source_arbitrary"
                        
            else:

                #################################################################
                ## Here, each different structure comes from one database so
                ## we have to pick one arbitrarily
                #################################################################

                #################################################################
                ## First, if it's a SMILE string and not InChI, it's difficult
                ## to parse, so just prioritize which database it came from
                ## and pick one
                #################################################################

                if(structure_type == "SMILE"):
                    if('MetaCyc' in sources_structures):
                        chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                        pick_reason="smile_db_priority:MetaCyc"
                    elif('KEGG' in sources_structures):
                        chosen_structure = sorted(sources_structures['KEGG'])[0]
                        pick_reason="smile_db_priority:KEGG"
                    elif('ChEBI' in sources_structures):
                        chosen_structure = sorted(sources_structures['ChEBI'])[0]
                        pick_reason="smile_db_priority:ChEBI"
                    elif('Rhea' in sources_structures):
                        chosen_structure = sorted(sources_structures['Rhea'])[0]
                        pick_reason="smile_db_priority:Rhea"

                else:

                    #################################################################
                    ## Secondly, if it's not a SMILE, then as InChI, we can break
                    ## the structure down and compare its connectivity which is
                    ## the first InChIKey 'layer'. But, we need to re-configure
                    ## this to use InChI
                    #################################################################

                    connected_structures = dict()
                    for structure in struct_conflicts:
                        connectivity = structure.split('-')[0]
                        if(connectivity not in connected_structures):
                            connected_structures[connectivity] = dict()
                        connected_structures[connectivity][structure]=1
                        
                    #################################################################
                    ## Having split the InChiKey into its connectivity layers and
                    ## grouping them, we take the connectivity that has more than
                    ## one structure (i.e. identical connectivity) but we're assuming
                    ## that there are not more than one different connectivities
                    ## that have multiple structures
                    #################################################################

                    chosen_connectivity = None
                    for connectivity in connected_structures:
                        if(len(connected_structures[connectivity])>1):
                            #At time of writing, only happens once per compound
                            chosen_connectivity = connectivity

                    #################################################################
                    ## Here, we have a connectivity that is the same for multiple
                    ## structures, so we take those structures and try to pick one
                    ## to use as the primary structure. First we try to see if there 
                    ## is one (and one only) that has stereochemistry. If not,
                    ## again, we prioritize by database
                    #################################################################

                    if(chosen_connectivity is not None):
                        
                        stereo_structures = dict()
                        for structure in connected_structures[chosen_connectivity]:
                            if('UHFFFAOYSA' not in structure):
                                stereo_structures[structure]=1

                        if(len(stereo_structures)==1):
                            chosen_structure=list(stereo_structures.keys())[0]
                            pick_reason="connectivity_stereochemistry"
                        else:
                            if('MetaCyc' in sources_structures):
                                chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                                pick_reason="connectivity_db_priority:MetaCyc"
                            elif('KEGG' in sources_structures):
                                chosen_structure = sorted(sources_structures['KEGG'])[0]
                                pick_reason="connectivity_db_priority:KEGG"
                            elif('ChEBI' in sources_structures):
                                chosen_structure = sorted(sources_structures['ChEBI'])[0]
                                pick_reason="connectivity_db_priority:ChEBI"
                            elif('Rhea' in sources_structures):
                                chosen_structure = sorted(sources_structures['Rhea'])[0]
                                pick_reason="connectivity_db_priority:Rhea"

                    #################################################################
                    ## Finally if we get to the point where there are multiple
                    ## connectivities that only have one structure, we prioritize
                    ## by database. There's some redundancy here that can be
                    ## eliminated
                    #################################################################

                    if(chosen_structure is None):
                        #Here we have structures with the same formula, but different connectivity
                        #For now, we pick MetaCyc
                        if('MetaCyc' in sources_structures):
                            chosen_structure = sorted(sources_structures['MetaCyc'])[0]
                            pick_reason="multi_connectivity_db_priority:MetaCyc"
                        elif('KEGG' in sources_structures):
                            chosen_structure = sorted(sources_structures['KEGG'])[0]
                            pick_reason="multi_connectivity_db_priority:KEGG"
                        elif('ChEBI' in sources_structures):
                            chosen_structure = sorted(sources_structures['ChEBI'])[0]
                            pick_reason="multi_connectivity_db_priority:ChEBI"
                        elif('Rhea' in sources_structures):
                            chosen_structure = sorted(sources_structures['Rhea'])[0]
                            pick_reason="multi_connectivity_db_priority:Rhea"

            #################################################################
            ## We have a chosen structure and we collect the aliases for that
            ## structure, so that we can then collect all the structure of
            ## different types for the same aliases
            #################################################################

            chosen_aliases=dict()
            for source in struct_conflicts[chosen_structure]:
                for alias in struct_conflicts[chosen_structure][source]:
                    chosen_aliases[alias]=1

            #################################################################
            ## Having collected aliases for the same structure we iterate
            ## through the types, and the aliases, to find the right structure
            ## to print with the alias
            #################################################################

            for structure_type in "SMILE","InChIKey","InChI":
                if(structure_type not in Structs):
                    continue

                structure_to_use = None
                for alias in chosen_aliases:
                    for structure in Structs[structure_type][struct_stage]:
                        if(alias in Structs[structure_type][struct_stage][structure]):
                            structure_to_use=structure

                #################################################################
                ## Now we collect *all* aliases that are associated with different
                ## structures for the same compound, to associate with the
                ## primary structure
                #################################################################

                aliases=dict()
                for structure in Structs[structure_type][struct_stage]:
                    for alias in Structs[structure_type][struct_stage][structure]:
                        aliases[alias]=1

                #################################################################
                ## Finally, write to file (SMILES canonicalized to keep Unique
                ## idempotent with Recanonicalize_SMILES.py)
                #################################################################

                structure_to_use_out = (CANONICALIZE_SMILES(structure_to_use)
                                        if structure_type == "SMILE" else structure_to_use)
                unique_structs_file.write("\t".join((msid,\
                                                     structure_type,\
                                                     ";".join(sorted(aliases)),\
                                                     formula_charge_dict['formula'],\
                                                     formula_charge_dict['charge'],\
                                                     structure_to_use_out))+"\n")
                                        
    #################################################################
    ## Here we report the structural conflicts
    #################################################################

    if(struct_conflict==1):
        for structure in Structs[struct_type][struct_stage]:
            for external_id in Structs[struct_type][struct_stage][structure]:
                structure_conflicts_file.write("\t".join((msid,struct_type,struct_stage,structure,external_id,
                                                          Structs[struct_type][struct_stage][structure][external_id]))+"\n")
                                        
    #################################################################
    ## Here we report the formula conflicts
    #################################################################

    if(formula_conflict==1):
        for formula in Formulas[struct_type][struct_stage]:
            for external_id in Formulas[struct_type][struct_stage][formula]:
                formula_dict=json.loads(formula)
                formula_conflicts_file.write("\t".join((msid,struct_type,struct_stage,formula_dict['formula'],formula_dict['charge'],external_id,
                                                        Formulas[struct_type][struct_stage][formula][external_id]))+"\n")

    #################################################################
    ## Record which branch of the picker fired for this compound.
    ## "single_structure" means no choice was needed; the *_db_priority
    ## reasons are fallbacks when no source agreement existed; the
    ## stereochemistry reasons are picks made to prefer stereo-specified
    ## structures over UHFFFAOYSA.  formula_conflict means no pick was
    ## made (the picker ran but rival structures had different formulas).
    #################################################################

    if(struct_pass and pick_reason is not None):
        chosen_aliases_list=[]
        try:
            chosen_aliases_iter = chosen_aliases  # set when conflict cascade ran
        except NameError:
            chosen_aliases_iter = None
        if(struct_conflict==0):
            # Single-structure case: aliases are all those for the picked struct
            for struct_t in ("InChI","InChIKey","SMILE"):
                if(struct_t in Structs and struct_stage in Structs[struct_t]):
                    for s in Structs[struct_t][struct_stage]:
                        for a in Structs[struct_t][struct_stage][s]:
                            chosen_aliases_list.append(a)
                    break
        else:
            if(chosen_aliases_iter is not None):
                chosen_aliases_list = list(chosen_aliases_iter.keys()) if isinstance(chosen_aliases_iter, dict) else list(chosen_aliases_iter)
        chosen_struct_str = ""
        try:
            if(chosen_structure is not None):
                chosen_struct_str = chosen_structure
        except NameError:
            pass
        if(not chosen_struct_str and struct_conflict==0):
            chosen_struct_str = list(Structs[struct_type][struct_stage].keys())[0]
        pick_reasons_file.write("\t".join((msid,struct_type,struct_stage,pick_reason,
                                           chosen_struct_str,
                                           ";".join(sorted(set(chosen_aliases_list)))))+"\n")
    elif(formula_conflict==1):
        pick_reasons_file.write("\t".join((msid,struct_type,struct_stage,"formula_conflict_no_pick","",""))+"\n")
