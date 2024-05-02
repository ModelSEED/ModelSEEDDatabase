# ModelSEED Biochemistry Database Scripts (Structures)

There are several scripts that we use to organize and collate the structures that we wish to use as part of the biochemistry:

* `Print_Structure_Formula_Charge.py`

This script takes all the InChI and SMILE strings from the original
sources and extracts the formulas and charges, including protonation
states You should only have to do this whenever you update the
structures themselves, but, you might get changes occuring if you
install updated versions of RDKit and/or OpenBabel in your conda
environment, so double-check. At time of submission we used RDKit
2020.03.1.0 and OpenBabel 2.4.1. I have an outline of what I need to 
do to get these working on a new mac below.

EDIT: As of 09/15/22 we used RDKit 2022.03.5 and OpenBabel 3.1.1
EDIT: This the same as of 01/23/24

* `List_ModelSEED_Structures.py` 

This script takes the list of InChI and SMILE strings from the
original sources and attempts to consolidate them for the ModelSEED
database.  As it checks for formula conflicts, it needs the output of
the previous script and it reports structure and formula conflicts in
files in the same directory.  Its two key output files are
`../../Biochemistry/Structures/All_ModelSEED_Structures.txt` and
`../../Biochemistry/Structures/Unique_ModelSEED_Structures.txt`.

The latter file is the structural heart, and our goal is to expand on
it via curation of the conflicts and integration of more structures.

* `Update_Compound_Structures_Formulas_Charge.py`

This script takes the output of the previous two scripts, and uses
them to update the ModelSEED database.

* `Update_Compound_pKas.py`

This script takes the pKa files found in
`../../Biochemistry/Structures`, and uses them to update the ModelSEED
database.

* `Report_Redundant_Structures.py`
* `Report_Conflicting_Structures.py`

These two scripts build tab-separated files of conflicts and
redundancies for reporting. Their output isn't necessarily optimal, or
include all possible problems, but is used as a starting point for
loading into spreadsheets with the goal of curating by a team.

Redundant structures are structures that are found assigned to more
than one compound, and could be merged. Concflicting structures are
multiple but different structures found assigned to a single compound
and could be disambiguated.

# Addendum

Acyl-carrier proteins are frequently found in metabolic
reconstructions and are a key component of some pathways such as fatty
acid biosynthesis. As the protein chains vary from species to species
there's no determined structure to use, and frequently the formula and
the charge of the phosphopantetheine prosthetic group as well as the
attached fatty acyl chain can be overlooked and leads to reaction
imbalance. Here we attempt to manually curate the formula and charge
that would maintain the mass-balance of the fatty acid biosynthetic
pathways, and others, in the `ACPs_Master_Formula_Charge.txt` file.

# Installing RDKit and OpenBabel

You might need swig and its python bindings but the conda packages may 
install fine:

`port install swig`
`port install swig-python`

You'll want to install anaconda and follow this process for creating
a specific environment for the chemoinformatics packages. This was
the approach with the least pain:

`conda create -c conda-forge -n msd-env rdkit`
`conda activate msd-env`
`conda install openbabel -c conda-forge`
