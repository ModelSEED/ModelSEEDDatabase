# Biochemistry
The purpose of this folder is to version control updates, changes, and corrections to the current Biochemistry object using the sources from ModelSEED, PlantSEED, and Probabilistic Annotation.  May there be one version to rule them all.

A Biochemistry object is defined by its compounds and reactions.  The following files contain the data needed to build a Biochemistry object.

* compounds.master.tsv: List of compounds for combined master biochemistry
* compounds.master.mods: Modifications to apply to combined master biochemistry
* compounds.default.tsv: List of compounds from ModelSEED biochemistry
* compounds.plantdefault.tsv: List of compounds from PlantSEED
* reactions.master.tsv: List of reactions for combined master biochemistry
* reactions.master.mods: Modifications to apply to combined master biochemistry
* reactions.default.tsv: List of reactions from ModelSEED biochemistry
* reactions.plantdefault.tsv: List of reactions from PlantSEED biochemistry

To build a Biochemistry object from the files, change to the scripts directory and run the following commands:

1. `./Print_Master_Compounds_List.pl` merges the default and plantdefault files, applies modifications, and creates the master file.
2. `./Print_Master_Reactions_From_Files.pl` merges the default and plantdefault files, applies modifications, checks reaction mass and charge balance, and creates the master file.
3. `./Build_Biochem_JSON.pl master-2015a` creates a Biochemistry object with the ID master-2015a and exports it to a JSON file.

## Compound file format
A compound file describes the compounds (or metabolites) involved in biochemical reactions.  There is one compound per line with fields separated by tabs.  The following fields are required.

1. **id**: Unique ID for the compound in the format cpdNNNNN where NNNNN is a five digit number (e.g. cpd00001)
2. **abbreviation**: Short name of compound
3. **name**: Long name of compound
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
17. **abstract_compound**: _Need definition_ or "null" if not specified (currently all compounds are set to null)
18. **comprised_of**: _Need definition_ or "null" if not specified (currently all compounds are set to null)
19. **aliases**: List of alternative names of compound separated by semicolon or "null" if not specified (see below for description of format)

### Format of pka and pkb
The pka and pkb fields are in this format:

    atoms:value

where "atoms" is the number of atoms and "value" is the dissociation constant value.  Multiple pkas or pkbs are separated by a semicolon.  For example, this is the pka for NAD:

    17:1.8;18:2.56;6:12.32;25:11.56;35:13.12

### Format of aliases
An alias is in this format:

	"source:value"
	
where "source" is the name of the alternative database and "value" is the name or ID in the alternative database. Multiple aliases are separated by a semicolon.  For example, this is the list of aliases for Cobamide (or cpd00181):

	"KEGG:C00210";"name:Cobamide";"searchname:cobamide";"ModelSEED:cpd00181";"KBase:kb|cpd.181"

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
3. **name**: Long name of reaction
4. **code**: Definition of reaction expressed using compound IDs and before protonation (see below for description of format)
5. **stoichiometry**: Definition of reaction expressed in stoichiometry format (see below for description of format)
6. **is_transport**: True if reaction is a transport reaction 
7. **equation**: Definition of reaction expressed using compound IDs and after protonation (see below for description of format)
8. **definition**: Definition of reaction expressed using compound names (see below for description of format)
9. **reversibility**: Reversibility of reaction where ">" means right directional, "<" means left directional, "=" means bi-directional, and "?" means unknown (_Need a better description of reversibility and direction_)
10. **direction**: Direction of reaction where ">" means right directional, "<" means left directional, and "=" means bi-directional
11. **abstract_reaction**: _Need definition_ or "null" if not specified (currently all reactions are set to null)
12. **pathways**: Pathways reaction is a part of or "null" if not specified (currently all reactions are set to null)
13. **aliases**: List of alternative names of reaction separated by semicolon or "null" if not specified (format is the same as Compounds file)
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

where "n" is the compound coefficient, "cpdid" is the compound ID, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]

### Format of reaction definition using compound names 
Each compound participating in the reaction is in this format:

	(n) cpdname[m]

where "n" is the compound coefficient, "cpdname" is the compound name, and "m" is a compartment index number.  Compounds are separated by "+" and reactant and product compounds are delimited by the direction symbol.  For example, this is the definition of reaction rxn00001:

	(1) H2O[0] + (1) PPi[0] <=>  (2) Phosphate[0] + (1) H+[0]

### Format of reaction stoichiometry
Each compound participating in the reaction is in this format:

	n:cpdid:m:i:"cpdname"

where "n" is the compound coefficient and a negative number indicates a reactant and a positive number indicates a product, "cpdid" is the compound ID, "m" is the compartment index number, "i" is the community index number (I think this is no longer needed), and "cpdname" is the compound name.  Compounds are separated by semicolon.  For example, this is the stoichiometry of reaction rxn00001:

	-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"

### Reaction status values

The reaction status field is a string with one or more values. Multiple values are separated with a "|" character.  There are multiple values when a reaction is updated after mass and charge balancing.  If "OK" is one of the values, the reaction is valid and additional values describe the changes made to the reaction to balance it.  The status values are:

* OK means the reaction is valid.  If "OK" is the only value, then the reaction was valid with no changes.
* MI means there is a mass imbalance. The remainder of the string after the first colon indicates what atoms are unbalanced and the number of atoms needed to balance the reaction.  Multiple atoms are separated by a "/" character.  A positive number means the righthand side of the reaction has more atoms and a negative number means the lefthand side of the reaction has more atoms.
* CI means there is a charge imbalance.  A positive number after the first colon means the righthand side of the reaction has a larger charge and a negative number means the lefthand side of the reaction has a larger charge.
* HB means the reaction has been balanced by adding hydrogen to it.
* EMPTY means reactants cancel out completely.
* CPDFORMERROR means at least one compound either has no formula or has an invalid formula.

For example, rxn00277 has this definition:

	(1) Glycine[0] <=> (1) HCN[0] or
	(1) C2H5NO2 <=> (1) CHN
	
and its status shows that it has a mass imbalance with 1 extra carbon, 4 extra hydrogen, 2 extra oxygen atoms on the lefthand side of the reaction:

	MI:C:-1/H:-4/O:-2

And rxn00008 has this definition:

	(2) H2O[0] <=> (1) H2O2[0] + (2) H+[0]

and its status shows it has a charge imbalance with the righthand side of the reaction having a larger charge:

	CI:2

## Reaction modifications file format
A reaction modification file describes modifications to make to the master reaction file.  There is no header line in the file. Each line has these fields:

1. **id**: ID of reaction to modify
2. **source**: Source of modification where "plantdefault" is for PlantSEED and "curated" is for manual curation
3. **action**: (1) Name of field to modify, (2) "replace" to replace a compound ID with another compound ID in all fields, (3) "priority" to select the reaction from a specific database, (4) "add" to add a compound to the reaction, (5) "remove" to remove the reaction from the master file, or (6) "coefficient" to update the coefficient of a compound 
4. **value**: (1) When action is name of field, the modified value of the specified field, (2) when action is "replace", original compound ID and new compound ID, (3) when action is "priority", the name of database to select the reaction from, (4) when action is "add", the compound ID, coefficient (negative number for reactant, positive number for product), and compartment index number, (5) when action is "remove", the value is ignored but must be specified, or (6) when action is "coefficient", the compound ID and new coefficient value

Note that a reaction ID can be repeated if there are multiple fields that need to be modified.

## Archived files

The following files were used to merge the ModelSEED and PlantSEED biochemistry and are saved for reference.

* compounds.plantdefault_obs.tsv: List of compounds from PlantSEED biochemistry (includes obsolete IDs)
* compounds.tsv: First attempt at merged compounds file, currently empty
* reactions.default.cf.tsv: List of reactions from ModelSEED biochemistry in compartment-free notation
* reactions.plantdefault.cf.tsv: List of reactions from PlantSEED biochemistry in compartment-free notation
* reactions.plantdefault_obs.tsv: List of reactions from PlantSEED biochemistry (includes obsolete IDs)
* reactions.plantdefault_obs.cf.tsv: List of reactions from PlantSEED biochemistry in compartment-free notation (includes obsolete IDs)
* reactions.tsv: First attempt at merged reactions file, currently empty
* compartments.plantdefault_obs.tsv: List of compartments from PlantSEED biochemistry (includes obsolete IDs and is not used)
* compartments.master.tsv: List of compartments for combined master biochemistry
* compartments.default.tsv: List of compartments from ModelSEED biochemistry
* compartments.plantdefault.tsv: Set of compartments from PlantSEED biochemistry
* Workspaces/KBaseTemplateModels.rxn: List of reactions and the model templates that include the reaction

