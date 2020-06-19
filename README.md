<!--[![Build
Status](https://travis-ci.org/ModelSEED/ModelSEEDDatabase.svg?branch=master)](https://travis-ci.org/ModelSEED/ModelSEEDDatabase)-->
# ModelSEED Biochemistry Database

<a href="https://www.biorxiv.org/content/10.1101/2020.03.31.018663v2">The
ModelSEED Database for the integration of metabolic annotations and
the reconstruction, comparison, and analysis of metabolic models for
plants, fungi, and microbes.</a>

## Abstract

For over ten years, ModelSEED has been a primary resource for the
construction of draft genome-scale metabolic models based on annotated
microbial or plant genomes. Now being released, the biochemistry
database serves as the foundation of biochemical data underlying
ModelSEED and KBase. The biochemistry database embodies several
properties that, taken together, distinguish it from other published
biochemistry resources by: (i) including compartmentalization,
transport reactions, charged molecules and proton balancing on
reactions;; (ii) being extensible by the user community, with all data
stored in GitHub; and (iii) design as a biochemical “Rosetta Stone” to
facilitate comparison and integration of annotations from many
different tools and databases. The database was constructed by
combining chemical data from many resources, applying standard
transformations, identifying redundancies, and computing thermodynamic
properties. The ModelSEED biochemistry is continually tested using
flux balance analysis to ensure the biochemical network is
modeling-ready and capable of simulating diverse
phenotypes. Ontologies can be designed to aid in comparing and
reconciling metabolic reconstructions that differ in how they
represent various metabolic pathways. ModelSEED now includes 33,978
compounds and 36,645 reactions, available as a set of extensible files
on GitHub, and available to search at https://modelseed.org and KBase.

## Related Links

* <a href="https://modelseed.org">ModelSEED</a>
* <a href="https://github.com/ModelSEED/ModelSEEDDatabase">(this) Git repository</a>
* <a href="https://kbase.us">KBase</a>

## Documentation

We provide some documentation for the database and everything required
in this repository to build and maintain it in these directories:

### [Biochemistry](Biochemistry/README.md)

This holds the main database files.

### [Scripts](Scripts/README.md)

This holds the scripts that we use to build and maintain the database.

### [Solr](Solr/README.md)

This holds the Solr schema and examples on how to access the
biochemistry at various endpoints.
