# ModelSEED Biochemistry Database Scripts (Structures)

There are several scripts that we use to organize and collate the structures that we wish to use as part of the biochemistry:

* `Print_Structure_Formula_Charge.py`

This script takes all the InChI and SMILE strings from the original
sources and extracts the formulas and charges, including protonation
states You should only have to do this whenever you update the
structures themselves, but, you might get changes occuring if you
install updated versions of RDKit and/or OpenBabel in your conda
environment, so double-check. At time of submission we used RDKit
2020.03.1.0 and OpenBabel 2.4.1

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

* `Update_Compound_Structures_Formulas_Charge.py`

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
