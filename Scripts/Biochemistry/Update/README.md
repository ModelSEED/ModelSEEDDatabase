# There's three groups of scripts here for updating the database

## A: Adding new biochemistry altogether
# documentation explaining the files needed will be forthcoming
./Add_New_Compounds.py
./Merge_Formulas.py
./Add_New_Reactions.py

## B: Adding or removing an attribute for a biochemical entity
# this will rarely be needed, but a user can extend or delimit the
# fields they want for their own uses
./Add_New_Biochemistry_Attribute.py
./Remove_Biochemistry_Attribute.py

## C: Removing obsolete entities after curation. When new biochemistry
# is added, but it is found that it is obsolete because it matches
# what is already in the database, then these scripts will make sure
# they are deleted
./Remove_Newly_Obsolescent_Compounds.py
./Remove_Newly_Obsolescent_Reactions.py

## NB: If any of these are run, all scripts in the Refresh folder must be run in order to make sure
## the entire database is up to date.