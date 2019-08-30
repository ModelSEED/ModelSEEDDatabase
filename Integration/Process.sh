#!/bin/bash
export PYTHONPATH=/Users/seaver/Projects/ModelSEEDDatabase/Libs/Python/
#git checkout ../Biochemistry/compounds.* ../Biochemistry/Aliases/Unique_ModelSEED_Compound_*.txt
../Scripts/Biochemistry/Add_New_Compounds.py Yeast_8.3.4_2019/yeastGEM_Reactions.tbl KEGG 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py Yeast_8.3.4_2019/yeastGEM_Reactions.tbl 'Published Model' f
../Scripts/Biochemistry/Add_New_Compounds.py iCY1106_25582171_2015/iCY1106_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iCY1106_25582171_2015/iCY1106_Reactions.tbl 'Published Model' f