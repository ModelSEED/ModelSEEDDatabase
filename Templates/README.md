# Model Templates
A Model Template is the source that is used to reconstruct a Model for an organism.  A Model Template describes the biomass and reactions for a type of organism.  A Model Template can be generic for a large set of organisms or specific to a small set of organisms. These are the current Model Templates:

* **GramNegative** for models of gram negative bacteria
* **GramPositive** for models of gram positive bacteria
* **Core** for models of core metabolism of bacteria
* **Plant** for models of plants

There is a folder for each Model Template which contains the following files:

* Compartments.tsv: List of compartments
* Biomasses.tsv: List of biomass reactions
* BiomassCompounds.tsv: List of compounds in the biomass reactions
* Reactions.tsv: List of reactions (subset of master biochemistry)

When a Model is reconstructed from the Model Template only the compartments, biomass reactions, and reactions from the above files are available to be included in the Model.

## Compartments file format
The Compartments file describes the compartments that are available for Models reconstructed from the Model Template.  There is one compartment per line with fields separated by tabs.  Each line has the following fields:

* **index**: Compartment index number from Biochemistry reaction stoichiometry
* **id**: Unique ID of compartment (not sure how to handle community models where compartments are numbered)
* **name**: Long name of compartment
* **hierarchy**: Number describing where compartment is in the hierarchy, starting from 0 for extracellular (is this used anywhere?)
* **pH**: Value of pH of compartment (is this used anywhere?)
* **aliases**: List of alternative names of compartment separated by semicolon or "null" if not specified

## Biomasses file format
The Biomasses file describes the biomass reactions. There is one biomass reaction per line with fields separated by tabs. The compounds involved in the biomass reaction are described in the Biomass Compounds file. Each line has the following fields:

* **id**: Unique ID of biomass reaction
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
The Biomass Compounds file describes the compounds that are included in the biomass reactions. There is one compound per line with fields separated by tabs. Each line has the following fields:

* **biomass_id**: ID of biomass reaction (must match an ID in the Biomasses file)
* **id**: ID of compound (must match an ID in master biochemistry)
* **coefficient**: Coefficient value where a negative number indicates a reactant and a positive value indicates a product
* **coefficient_type**: Type of coefficient (valid values are "MOLFRACTION", "MOLSPLIT", "EXACT", "AT", "GC", or "MULTIPLIER")
* **class**: Class of compound (valid values are "dna", "rna", "cellwall", "lipid", "protein", "cofactor", "energy", or "other")
* **linked_compounds**: List of linked compounds and coefficients or "null" if not specified (see below for description of format) 
* **compartment**: ID of compartment where biomass reaction occurs (Could it ever be a compartment other than cytosol?)

### Linked compounds format
A linked compound is in this format:

	cpdid:coeff

where "cpdid" is the ID of the linked compound and "coeff" is the coefficient for the compound.  A negative coefficient means X and a positive coefficient means Y.  Multiple compounds are separated by a "|" character.

## Reactions file format
The Reactions file describes the reactions that can be added to a Model reconstructed from a Model Template. There is one reaction per line with fields separated by tabs.  Each line has the following fields:

* **id**: ID of reaction (must match an ID in master biochemistry)
* **compartment**: List of compartment IDs separated by "|" (see below for detailed description)
* **direction**: Default direction of the reaction when added to a Model by gene association (valid values are "<", "=", ">")
* **gfdir**: Direction of the reaction when directionality is reversed by gap fill (can be unspecified) (not seeing any reactions with field set)
* **type**: Type of the reaction (see below for detailed description)
* **base_cost**: Base cost to add reaction to a Model (encodes all penalities)
* **forward_cost**: Cost to add reaction in forward direction to a Model
* **reverse_cost**: Cost to add reaction in reverse direction to a Model
* **complexes**: List of complex IDs separated by "|" character or null if not specified

### Compartment field description
Each compartment ID in the list must match the ID of a compartment from the compartments file.  The compartment IDs are used to replace the compartment index in the definition of the reaction from the master biochemistry. The first compartment ID in the list is substituted for compartment index 0, the second compartment ID is substituted for compartment index 1, and so forth.

In addition, the first compartment ID in the list is used as the suffix for the template reaction ID.

### Type field description
The type field has the following valid values:

* **conditional**: A conditional reaction can be added to a Model only if the associated gene is in the organism
* **gapfilling**: A gapfilling reaction can be added to a Model as needed (could go away if we get gene associations)
* **spontaneous**: A spontaneous reaction can be added to a Model even if the associated gene is not in the organism
* **universal**: An universal reaction is always added to a Model

## Complexes file format
cpxid
cpxname
reference Where complex came from (PUBMED id of publication) (could get from KBase DLITs for SEED?)
source KEGG published or SEED
confidence Value from 0 to 1 
optional Is the role necessary?
triggering If role is there does it trigger reaction to be there
type
roleid
rolename
roletype 
rolefeatures Reference to feature in blast database

The Complexes.tsv file contains