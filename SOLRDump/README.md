# Files for importing to SOLR
The purpose of this folder is to hold the files used for importing the databases into SOLR.  There are these files:

* Compounds.tsv: Master compounds file
* Reactions.tsv: Master reactions file
* TemplateReactions.tsv: Template reactions file

## Compounds file format
The Compounds file describes the compounds from the master biochemistry. There is one compound per line with fields separated by tabs. Each line has the same fields as the master compounds file.

## Reactions file format
The Reactions file describes the reactions from the master biochemistry. There is one reaction per line with fields separated by tabs. Each line has the same fields as the master reactions file.

## Template Reactions file format
The Template Reactions file describes the reactions that are included in Model Templates.  There is one reaction per line with fields separated by tabs.  Each line has the following fields:

* **id**: Unique ID three parts, template id, reaction id, compartment id
* **reaction_id**: ID for the reaction in the master biochemistry in the format rxnNNNNN where NNNNN is a five digit number (e.g. rxn03789)
* **abbreviation**: Short name of reaction
* **name**: Human readable long name of reaction
* **code**: Definition of reaction expressed using compound IDs and before protonation (see below for description of format)
* **stoichiometry**: Definition of reaction expressed in stoichiometry format (see below for description of format)
* **is_transport**: True if reaction is a transport reaction
* **equation**: Definition of reaction expressed using compound IDs and after protonation (see below for description of format)
* **definition**: Definition of reaction expressed using compound names (see below for description of format)
* **model_direction**: Default direction when reaction is added to model by gene association
* **gapfill_direction**: Direction when gap fill reverses directionality
* **type**: What?  Seems to be null 
* **base_cost**: Cost to add reaction to a model (encodes penalties)
* **forward_penalty**: Cost to add reaction to a model in forward direction
* **reverse_penalty**: Cost to add reaction to a model in reverse direction
* **pathways**: List of pathways
* **aliases**: List of alternative names of compound separated by semicolon or "null" if not specified
* **ec_numbers**: List of EC numbers
* **deltag**: Value for change in free energy of compound or "null" when unknown
* **deltagerr**: Value for change in free energy error of compound or "null" when unknown
* **template_id**: Unique ID of model template
* **template_name**: Name of model template
* **template_modeltype**: Type of model template (valid values are "genome_scale_model" or "core_model")
* **template_domain**: Domain of organisms (valid values are "Bacteria" or "Plant") 
* **template_version**: Version number of model template (not sure how this is ever set in model template)
* **template_is_current**: True if model template is current (not sure how this is determined)
* **template_owner**: User name of template owner (not sure this is really needed -- seems to be specific to workspace objects)
* **compartment_ids**: List of compartment IDs in format described below (0:c)
* **complex_ids**: List of complex IDs separated by semicolon for what?
* **compound_ids**: List of compound IDs separated by semicolon of compounds involved in biomass reaction

## Template Biomasses file
The Template Biomasses file describes the biomass that is included in Model Templates.  There is one biomass per line with fields separated by tabs.  Each line has the following fields:

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
* **template_id**: Unique ID of model template
* **template_name**: Name of model template
* **template_modeltype**: Type of model template (valid values are "genome_scale_model" or "core_model")
* **template_domain**: Domain of organisms (valid values are "Bacteria" or "Plant")
* **template_version**: Version number of model template (not sure how this is ever set in model template)
* **template_is_current**: True if model template is current (not sure how this is determined)
* **template_owner**: User name of template owner (not sure this is really needed -- seems to be specific to workspace objects)
* **compartment_ids**: List of compartment IDs in format described below (0:c)
* **compound_ids**: List of compound IDs separated by semicolon of compounds involved in biomass reaction
* **compound_data**: List of data separated by semicolon describing each compound's contribution to biomass reaction (see below for description of format)

### Compound data format

ID, name, coefficient, coefficient type, class, linked compounds, linked compound coefficient
cpd00002:"ATP":-0.262:MOLFRACTION:rna:cpd00012{-1}|

coefficient type: MOLFRACTION, MOLSPLIT, EXACT
