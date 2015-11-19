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

## Roles file format
The Roles file describes functional roles that are available for Models.  There is one master Roles file and a Model Template only includes the roles that are associated with the reactions in the Model Template.  There is one role per line with fields separated by tabs.  Each line has the following fields:

* **id**: ID of role
* **name**: Long name of role
* **source**: Source of role (valid values are "ModelSEED", "PlantSEED", "KEGG", or "published")
* **features**: List of features associated with role separated by semicolon or "null" if not specified
* **aliases**: List of aliases for role separated by semicolon or "null" if not specified

## Complexes file format
The Complexes file describes the complexes that are available for Models.  There is one master Complexes file and a Model Template only includes the complexes that are associated with the reactions in the Model Template.  There is one complex per line with fields separated by tabs.  Each line has the following fields:

* **id**: ID of complex
* **name**: Name of complex (typically the same as the ID)
* **source**: Source of complex (valid values are "ModelSEED", "PlantSEED", "KEGG", or "published")
* **reference**: Reference to where complex came from (_PUBMED id of publication or KBase DLITs for SEED?_) or "null" if not specified
* **confidence**: Value from 0 to 1 (_need details on how value is calculated_)
* **roles**: List of roles and relationship to the complex (format described below)

### Roles field format

* **role_id**: ID of role
* **type**: Type (always seems to be "role_mapping")
* **optional**: Is the role necessary?
* **triggering**: If role is found does it trigger associated reaction


## Compartments file format
The Compartments file describes the compartments that are available for Models reconstructed from the Model Template.  There is one compartment per line with fields separated by tabs.  Each line has the following fields:

* **index**: Compartment index number from Biochemistry reaction stoichiometry
* **id**: Unique ID of compartment (not sure how to handle community models where compartments are numbered)
* **name**: Long name of compartment
* **hierarchy**: Number describing where compartment is in the hierarchy, starting from 0 for extracellular (_Is hierarchy ever used in models?_)
* **pH**: Value of pH of compartment (_Is pH ever used in models?_)
* **aliases**: List of alternative names of compartment separated by semicolon or "null" if not specified (_Are aliases ever used in models?_)

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
* **coefficient_type**: Type of coefficient (valid values are "MOLFRACTION", "MOLSPLIT", "EXACT", "AT", "GC", or "MULTIPLIER") (need better description of these values)
* **class**: Class of compound (valid values are "dna", "rna", "cellwall", "lipid", "protein", "cofactor", "energy", or "other")
* **linked_compounds**: List of linked compounds and coefficients or "null" if not specified (see below for description of format) 
* **compartment**: ID of compartment where biomass reaction occurs (Could it ever be a compartment other than cytosol?)

### Linked compounds format
A linked compound is in this format:

	cpdid:coeff

where "cpdid" is the ID of the linked compound and "coeff" is the coefficient for the compound.  A negative coefficient means X and a positive coefficient means Y.  Multiple compounds are separated by a "|" character. _Need more information on how linked compounds are used in models_.

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

* **conditional**: A conditional reaction can be added to a Model only if the associated gene is in the organism. At least one complex ID is required for conditional reactions.
* **gapfilling**: A gapfilling reaction can be added to a Model as needed (could go away if we get gene associations).
* **spontaneous**: A spontaneous reaction can be added to a Model even if the associated gene is not in the organism. (is a complex required for spontaneous reactions?)
* **universal**: An universal reaction is always added to a Model.

## How to create a model template

#### Step 1
Create a folder for the template in the Templates folder and add a compartments file, biomasses file, biomass compounds file, and reactions file using the formats described above.

#### Step 2
Run Build\_Model\_Template.py to build a Model Template object from the source files and store it in a workspace.

## How to create complex and role files

Note that these steps only need to be run to generate the complex and role files from the original sources. In the future, the complex and roles can be updated directly (e.g. when KEGG adds new enzymes).

#### Step 1
Run Print\_ModelTemplate\_from\_Workspaces.pl script which creates a Mappings folder with subfolders for each of the Mapping objects in the current system. Each subfolder has Mapping\_Complexes.txt and Mapping\_Roles.txt files with the corresponding data from the Mapping object.  Note this step uses the KBase workspace and requires a KBase runtime environment.

#### Step 2
Run Print\_KEGG\_Complex\_Role.py script which creates complexes and roles files from a KEGG enzyme database.  You must have a source file in the KEGG file format. Note this step uses the uses packages from KBase and requires a KBase runtime environment.

#### Step 3
Run Build\_Complex\_File.py which merges complexes from all source files and creates a Complexes file as described above.

#### Step 4
Run Build\_Role\_File.py which merges all roles from all source files and creates a Roles file as described above.

