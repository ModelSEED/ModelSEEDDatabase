#!/bin/bash
export PYTHONPATH=/Users/seaver/Projects/ModelSEEDDatabase/Libs/Python/
git checkout ../Biochemistry/compounds.* ../Biochemistry/Aliases/Unique_ModelSEED_Compound_*.txt
#Yeast Consensus Model (8.3.4)
../Scripts/Biochemistry/Add_New_Compounds.py Yeast_8.3.4_2019/yeastGEM_Compounds.tbl KEGG 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py Yeast_8.3.4_2019/yeastGEM_Reactions.tbl 'Published Model' f
#iCY1106
../Scripts/Biochemistry/Add_New_Compounds.py iCY1106_25582171_2015/iCY1106_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iCY1106_25582171_2015/iCY1106_Reactions.tbl 'Published Model' f
#iJL1454
../Scripts/Biochemistry/Add_New_Compounds.py iJL1454_23624532_2013/iJL1454_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iJL1454_23624532_2013/iJL1454_Reactions.tbl 'Published Model' f
#iCT646
../Scripts/Biochemistry/Add_New_Compounds.py iCT646_26915092_2016/iCT646_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iCT646_26915092_2016/iCT646_Reactions.tbl 'Published Model' f
#iYL619
../Scripts/Biochemistry/Add_New_Compounds.py iYL619_PCP_23236514_2012/iYL619_PCP_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iYL619_PCP_23236514_2012/iYL619_PCP_Reactions.tbl 'Published Model' f
