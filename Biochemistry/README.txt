This document sits in the biochemistry folder of the ModelSEED/ModelSEEDDatabase repo on github. The hard link for this repo is https://github.com/ModelSEED/ModelSEEDDatabase.

PURPOSE: To version control updates and fixes to the latest biochemistry object for the ModelSEED, PlantSEED, and ProbAnno. May there be one version to rule them all.

ORGANIZATION: 

FOLDERS: 
(root) 
biochemistry -> folder containing a set of tab delimited tables with biochemistry data
(first children)
scripts -> folder with scripts for compiling tables into typed objects

FILES:
(root)
compounds.tsv -> description of the compounds data that describes the metabolites involved in biochemical reactions. 
COLUMNS: 
PrimaryName: Human readable compound name
Formula: Hill standard chemical formula in protonated form (to match the reported charge)
Charge: Electric charge of the compound
Structure (inchi/smiles):  ID corresponding with ID of structure in compounds.inchi.tsv
ID: 5 digit number associated with the compound – e.g., cmp00000
AliasIDs: semi-colon separated list of Ids for duplicate, outdated, or otherwise repeated compounds

compounds.inchi.tsv -> file of chemical structures that correspond to in compounds.tsv
COLUMNS:
ID: Unique ID assigned by ModelSEED
Inchi: Struture data in Inchi format
Source_id: ID assigned by source (KEGG, MetaCyc)
Source_db: Source (KEGG, MetaCyc)

reactions.tsv -> PENDING FUTURE DISCUSSION