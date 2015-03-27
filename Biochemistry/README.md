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

* id: Unique ID for the compound in the format cpdNNNNN where NNNNN is a five digit number (e.g. cpd00001)
* primary_name: Human readable compound name
* formula: Standard chemical format using Hill format in protonated form to match reported charge
* charge: Electric charge of compound

The following fields are optional:

* aliases: List of IDs for duplicate, out-dated or otherwise repeated compounds

## Reaction file format
A reaction file describes the biochemical reactions.  There is one reaction per line with fields separated by tabs.  The following fields are required:

* id: Unique ID for the reaction in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
* primary_name: Human readable reaction name
* equation: Equation expressed in compound IDs
* definition: Equation expressed in compound names
* status: String describing status of the reaction where OK means valid, MI means there is a mass imbalance, CI means there is a charge imbalance, HI means there is a hydrogen imbalance.
