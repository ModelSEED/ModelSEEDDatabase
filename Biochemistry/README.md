# Biochemistry
The purpose of this folder is to version control updates, changes, and corrections to the current Biochemistry object for ModelSEED, PlantSEED, and Probabilistic Annotation.  May there be one version to rule them all.

A Biochemistry object is defined by its compounds and reactions.  The following files contain the data needed to build a Biochemistry object.

* compounds.default.tsv: Set of base compounds from ModelSEED
* compounds.plantdefault.tsv: Set of additional compounds from PlantSEED
* compounds.tsv: Merged compounds file
* reactions.default.tsv: Set of base reactions from ModelSEED
* reactions.plantdefault.tsv: Set of additional reactions from PlantSEED
* reactions.tsv: Merged reactions file

See the scripts folder for commands to compile the tables into typed objects.

## Compound file format
A compound file describes the compounds (or metabolites) involved in biochemical reactions.  There is one compound per line with fields separated by tabs.  The following fields are required.

* **id**: Unique ID for the compound in the format cpdNNNNN where NNNNN is a five digit number (e.g. cpd00001)
* **name**: Human readable long name of compound
* **abbreviation**: Short name of compound
* **formula**: Standard chemical format (using Hill system) in protonated form to match reported charge
* **charge**: Electric charge of compound
* **isCofactor**: True if the compound is a cofactor  This is false for all compounds

The following fields are optional.

* **aliases**: List of IDs for duplicate, out-dated or otherwise repeated compounds or ID of compound in other databases
* **deltaG**: Change in free energy
* **deltaG_err**: 
* **mass**: Mass of compound (units?)
* **pkas**: do not understand the definition of the array
* **pkbs**: do not understand the definition of the array

## Reaction file format
A reaction file describes the biochemical reactions.  There is one reaction per line with fields separated by tabs.  The following fields are required:

* **id**: Unique ID for the reaction in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
* **name**: Human readable long name of reaction
* **abbreviation**: Short name of reaction
* **direction**: Direction of reaction where ">" means right directional, "<" means left directional, and "=" means bi-directional
* **thermoReversibility**: Reversibility of reaction where ">" means right directional, "<" means left directional, "=" means bi-directional, and "?" means unknown
* **defaultProtons**: Number of protons ???  This is 0 in all reactions
* **status**: String describing status of the reaction with multiple values delimited with a "|" character.  See below for details.
* **equation**: Equation expressed in compound IDs

What about the following fields?

* aliases: List of IDs for duplicate, out-dated or otherwise repeated reactions or ID of reaction in other databases
* deltaG: Change in free energy
* deltaG_err:

### Reaction status values ###

* OK means the reaction was valid with no changes.
* MI means there was a mass imbalance that was corrected.
* CI means there was a charge imbalance that was corrected.
* HI means there was a hydrogen imbalance?

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
