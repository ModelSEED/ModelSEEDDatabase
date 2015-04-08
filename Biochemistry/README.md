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
* **isCofactor**: True if the compound is a cofactor

The following fields are optional.

* **aliases**: List of IDs for duplicate, out-dated or otherwise repeated compounds or ID of compound in other databases
* **deltaG**: Change in free energy
* **deltaG_err**: 
* **mass**: Mass of compound (units?)
* **pkas**: do not understand the definition of the array
* **pkbs**: do not understand the definition of the array

## Reaction file format
A reaction file describes the biochemical reactions.  There is one reaction per line with fields separated by tabs.  The following fields are required:

* id: Unique ID for the reaction in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
* primary_name: Human readable reaction name
* equation: Equation expressed in compound IDs
* definition: Equation expressed in compound names
* status: String describing status of the reaction where OK means valid, MI means there is a mass imbalance, CI means there is a charge imbalance, HI means there is a hydrogen imbalance. Multiple values are delimited with a | character.

What about the following fields?

* aliases: List of IDs for duplicate, out-dated or otherwise repeated reactions or ID of reaction in other databases
* abbreviation: Short name of reaction
* deltaG: Change in free energy
* deltaG_err:
* reversibility: Need definition
