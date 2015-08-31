# Model Templates

A Model Template is the source that is used as input to reconstruct a Model for an organism.  A Model Template describes the biomass and reactions for a type of organism.  A Model Template can be generic for a large set of organisms or specific to a small set of organisms. These are the current Model Templates:

* GramNegative for models of gram negative bacteria
* GramPositive for models of gram positive bacteria
* Core for models of core metabolism of bacteria
* Plant for models of plants

There is a folder for each Model Template which contains the following files:

* BiomassCompounds.tsv: Compounds in the biomass reaction
* Reactions.tsv: Reactions a Model reconstructed from the Model Template can only include these reactions (subset of Biochemistry)

## Biomasses file format
The Biomasses file describes a biomass reaction. There is one biomass per line with fields separated by tabs.  Each line has the following fields:

* **id**: Unique ID 3 parts
* **name**: Name of biomass
* **type**: Type of biomass (valid values are "growth")
* **other**: Amount of other components in biomass (units?)
* **dna**: Amount of DNA in biomass (units?)
* **rna**: Amount of RNA in biomass (units?)
* **protein**: Amount of protein in biomass (units?)
* **lipid**: Amount of lipids in biomoass (units?)
* **cellwall**: Amount of cellwall in biomass (units?)
* **cofactor**: Amount of cofactors in biomass (units?)
* **energy**: Amount of energy in biomass (units?)

## Biomass Compounds file format
The Biomass Compounds file describes the compounds that are included in the biomass reaction. There is one compound per line with fields separated by tabs. Each line has the following fields:

* **biomass_id**: ID of biomass reaction
* **id**: ID of compound
* **coefficient**: Coefficient value where a negative number indicates a reactant and a positive value indicates a product
* **coefficient_type**: Type of coefficient (valid values are "MOLFRACTION", "MOLSPLIT", "EXACT", "AT", "GC", or "MULTIPLIER")
* **class**: Class of compound (valid values are "dna", "rna", "cellwall", "lipid", "protein", "cofactor", "energy", or "other")
* **linked_compounds**: List of linked compounds and coefficients or "null" if not specified (see below for description of format) 
* **compartment**: ID of compartment ?? not sure on this one

### Linked compounds format
A linked compound 

	cpdid:coeff

where "cpdid" is the ID of the linked compound and "coeff" is the coefficient for the compound.

The Complexes.tsv file contains