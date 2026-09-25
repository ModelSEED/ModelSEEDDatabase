# ModelSEED Biochemistry Database Scripts (Biochemistry)

Here we keep the scripts for loading and manipulating the biochemistry
database, we separate them into several categories depending on how we
use them. The "parent" script, which does the simple job of loading
the data from the original `*.tsv` files and some supplemental files,
then re-printing them, is the one we re-run every time to double-check
that everything is in place:

### `./Reprint_Biochemistry.py`

(What we look for is if this script runs silently, and `git status -s` shows no change)
 
The categories in which we sort the other scripts are:

## Updating

* `Add_New_Biochemistry_Attribute.py`
* `Remove_Biochemistry_Attribute.py`
* `Add_New_Compounds.py`
* `Add_New_Reactions.py`
* `Remove_Newly_Obsolescent_Compounds.py`
* `Remove_Newly_Obsolescent_Reactions.py`

## Refreshing
* `Adjust_Reaction_Protons.py`
* `Adjust_Reaction_Water.py`
* `Merge_Formulas.py`
* `Merge_Obsolete_Aliases.py`
* `Merge_Reactions.py`
* `Rebalance_Reactions.py` (pass `save` to write the recomputed statuses
  back; without it the script is a dry run. Its argument guard rejected
  `save` until 2026-09-24, so the rebalance step of `Refresh_Reactions.sh`
  had not been running.)
* `Rebuild_Reactions.py`
* `Update_Compound_Aliases.py`
* `Update_Reaction_Aliases.py`
* `Remove_Duplicate_Aliases.py`

## Maintaining
* `Check_Charges.py`
* `Check_Formulas.py`
* `Check_Links.py`
* `Check_Template_Reactions.py`
* `Check_Transport.py`
* `Fix_Values.py`
* `Manual_Update_Links.py`