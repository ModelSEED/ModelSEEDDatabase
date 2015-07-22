# Biochemistry
The purpose of this folder is to version control updates, changes, and corrections to the current Biochemistry object for ModelSEED, PlantSEED, and Probabilistic Annotation.  May there be one version to rule them all.

A Biochemistry object is defined by its compounds, reactions, and compartments.  The following files contain the data needed to build a Biochemistry object.

* compounds.master.tsv: List of compounds for combined master biochemistry
* compounds.master.mods: Modifications to apply to combined master biochemistry
* compounds.default.tsv: List of compounds from ModelSEED biochemistry
* compounds.plantdefault.tsv: List of compounds from PlantSEED

* reactions.master.tsv: List of reactions for combined master biochemistry
* reactions.master.mods: Modifications to apply to combined master biochemistry
* reactions.default.tsv: List of reactions from ModelSEED biochemistry
* reactions.plantdefault.tsv: List of reactions from PlantSEED biochemistry
* Workspaces/KBaseTemplateModels.rxn: List of reactions and the model templates that include the reaction 

* compartments.master.tsv: List of compartments for combined master biochemistry
* compartments.default.tsv: List of compartments from ModelSEED biochemistry
* compartments.plantdefault.tsv: Set of compartments from PlantSEED biochemistry

Run the following commands to compile the tables into typed objects.

1. scripts/Print\_Master\_Compounds\_List.pl merges the default and plantdefault files, applies modifications, and creates the master file.
2. scripts/Print\_Master\_Reactions\_List.pl merges the default and plantdefault files, applies modifications, and creates the master file.
3. scripts/Update\_Reaction\_Status.pl checks reaction mass and charge balance and updates reactions used in model templates.
4. scripts/Build\_Biochem\_JSON.pl creates a Biochemistry object and exports to a JSON file.

## Compound file format
A compound file describes the compounds (or metabolites) involved in biochemical reactions.  There is one compound per line with fields separated by tabs.  The following fields are required.

1. **id**: Unique ID for the compound in the format cpdNNNNN where NNNNN is a five digit number (e.g. cpd00001)
2. **abbreviation**: Short name of compound
3. **name**: Human readable long name of compound
4. **formula**: Standard chemical format (using Hill system) in protonated form to match reported charge
5. **mass**: Mass of compound or "null" when unknown
6. **source**: Source database of compound (currently only source is ModelSEED)
7. **structure**: Structure of compound using IUPAC International Chemical Identifier (InChI) format
8. **charge**: Electric charge of compound
9. **is_core**: True if compound is in core biochemistry (currently all compounds are set to true)
10. **is_obsolete**: True if compound is obsolete and replaced by different compound (currently all compounds are set to false)
11. **linked_compound**: List of compound IDs separated by semicolon related to this compound or "null" if not specified (used to link an obsolete compound to replacement compound)
12. **is_cofactor**: True if compound is a cofactor (currently all compounds are set to false)
13. **deltag**: Value for change in free energy of compound or "null" when unknown
14. **deltagerr**: Value for change in free energy error of compound or "null" when unknown
15. **pka**: Acid dissociation constants of compound (see below for description of format)
16. **pkb**: Base dissociation constants of compound (see below for description of format)
17. **abstract_compound**: Not sure of definition or "null" if not specified (currently all compounds are set to null)
18. **comprised_of**: Not sure of definition or "null" if not specified (currently all compounds are set to null)
19. **aliases**: Alternative names of compound or "null" if not specified (currently all compounds are set to null)

### Format of pka and pkb

The pka and pkb fields are in this format:

    atoms:value

where "atoms" is the number of atoms and "value" is the dissociation constant value.  Multiple pkas or pkbs are separated by a semicolon.  For example, this is the pka for NAD:

    17:1.8;18:2.56;6:12.32;25:11.56;35:13.12

## Compound modifications file format
A compound modification file describes modifications to make to the master compound file.  There is no header line in the file.  Each line has these fields:

1. **id**: ID of compound to modify
2. **source**: Source of modification where "plantdefault" is for PlantSEED and "curated" is for manual curation
3. **field**: Name of field to modify
4. **value**: Modified value of the specified field

Note that a compound ID can be repeated if there are multiple fields that need to be modified.

## Reaction file format
A reaction file describes the biochemical reactions.  There is one reaction per line with fields separated by tabs.  The following fields are required:

1. **id**: Unique ID for the reaction in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
2. **abbreviation**: Short name of reaction
3. **name**: Human readable long name of reaction
4. **code**: Definition of reaction expressed using compound IDs and before protonation (see below for description of format)
5. **stoichiometry**: Definition of reaction expressed in stoichiometry format (see below for description of format)
6. **is_transport**: True if reaction is a transport reaction 
7. **equation**: Definition of reaction expressed using compound IDs and after protonation (see below for description of format)
8. **definition**: Definition of reaction expressed using compound names (see below for description of format)
9. **reversibility**: Reversibility of reaction where ">" means right directional, "<" means left directional, "=" means bi-directional, and "?" means unknown
10. **direction**: Direction of reaction where ">" means right directional, "<" means left directional, and "=" means bi-directional
11. **abstract_reaction**: Not sure of definition or "null" if not specified (currently all reactions are set to null)
12. **pathways**: Pathways reaction is a part of or "null" if not specified (currently all reactions are set to null)
13. **aliases**: Alternative names of reaction or "null" if not specified (currently all reactions are set to null)
14. **ec_numbers**: Enzyme Commission numbers of enzymes that catalyze reaction or "null" if not specified (currently all reactions are set to null)
15. **deltag**: Value for change in free energy of reaction or 10000000 when unknown
16. **deltagerr**: Value for change in free energy error of reaction or 10000000 when unknown
17. **compound_id**: List of compound IDs separated by semicolon for compounds involved in reaction
18. **status**: String describing status of the reaction with multiple values delimited with a "|" character.  See below for details.
19. **is_obsolete**: True if reaction is obsolete and replaced by different reaction
20. **linked_reaction**: List of reaction IDs separated by semicolon related to this reaction or "null" if not specified (used to link an obsolete reaction to replacement reaction)

### Format of reaction definition using compound IDs 
Each compound participating in the reaction is in this format:

	(n) cpdid[m]

where "n" is the number of compounds, "cpdid" is the compound ID, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]

### Format of reaction definition using compound names 
Each compound participating in the reaction is in this format:

	(n) cpdname[m]

where "n" is the number of compounds, "cpdname" is the compound name, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) H2O[0] + (1) PPi[0] <=>  (2) Phosphate[0] + (1) H+[0]

### Format of reaction stoichiometry
Each compound participating in the reaction is in this format:

	n:cpdid:c:i:"cpdname"

where "n" is the number of compounds and a negative number indicates a reactant and a positive number indicates a product, "cpdid" is the compound ID, "c" is the compartment, "i" is the compartment index, and "cpdname" is the compound name.  Compounds are separated by semicolon.  For example, this is the stoichiometry of reaction rxn00001:

	-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"

### Reaction status values

* OK means the reaction is valid.  If "OK" is the only value, then the reaction was valid with no changes.
* MI means a mass imbalance was corrected.
* CI means a charge imbalance was corrected.
* HI means a hydrogen imbalance was corrected?
* HB means ?
* RC means the reversibility was corrected.
* FO means ?
* RO means ?
* SP means spontaneous?
* UN means ?

## Reaction modifications file format
A reaction modification file describes modifications to make to the master reaction file.  There is no header line in the file. Each line has these fields:

1. **id**: ID of reaction to modify
2. **source**: Source of modification where "plantdefault" is for PlantSEED and "curated" is for manual curation
3. **action**: (1) Name of field to modify, (2) "replace" to replace a compound ID with another compound ID in all fields, or (3) "priority" to select the reaction from a specific database
4. **value**: (1) When action is name of field the modified value of the specified field, (2) when action is "replace" original compound ID and new compound ID, or (3) when action is "priority" the name of database to select the reaction from

Note that a reaction ID can be repeated if there are multiple fields that need to be modified.

## Compartment file format
A compartment file describes the compartments in cells where biochemical reactions take place.  There is one compartment per line with fields separated by tabs:  The following fields are required:

* **id**: Unique ID for the compartment (single character)
* **name**: Human readable long name of compartment
* **hierarchy**: Number describing where compartment is in the hierarchy, starting from 0 for extracellular

## Complex role file format
A complex role file maps complexes to functional roles.  There is one complex role mapping per line with fields separated by tabs.  The following fields are required:

* **complex_id**: Unique ID for the complex in the format cpx.N where N is a number
* **complex_name**: Human readable name of complex (currently same as the ID)
* **complex_source**: Source of complex, valid values are ModelSEED and KEGG
* **complex_type**: Type of complex, valid values are SEED\_role\_complex and KEGG\_role\_complex
* **role_id**: Unique ID for the role in the format in the format fr.N where N is a number
* **role_name**: Human readable name of the role
* **role_source**: Source of role, valid values are SEED and KEGG
* **role_type**: Type of role, valid values are SEED\_role and KEGG\_role
* **role_aliases**: List of aliases for the role where each entry in the format type:alias
* **role_exemplar**: What is this?  Currently not set for any mappings
* **type**: What is this? All entries are set to role_mapping
* **triggering**: What is this? All entries are set to true
* **optional**: What is this? 12 of 2221 are set to true

## Obsolete files

The following files are obsolete and are saved for reference.

* compounds.plantdefault_obs.tsv: List of compounds from PlantSEED biochemistry (includes obsolete IDs)
* compounds.tsv: First attempt at merged compounds file, currently empty
* reactions.default.cf.tsv: List of reactions from ModelSEED biochemistry in compartment-free notation
* reactions.plantdefault.cf.tsv: List of reactions from PlantSEED biochemistry in compartment-free notation
* reactions.plantdefault_obs.tsv: List of reactions from PlantSEED biochemistry (includes obsolete IDs)
* reactions.plantdefault_obs.cf.tsv: List of reactions from PlantSEED biochemistry in compartment-free notation (includes obsolete IDs)
* reactions.tsv: First attempt at merged reactions file, currently empty
* compartments.plantdefault_obs.tsv: List of compartments from PlantSEED biochemistry (includes obsolete IDs and is not used)

