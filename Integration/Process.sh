#!/bin/bash
export PYTHONPATH=/Users/seaver/Projects/ModelSEEDDatabase/Libs/Python/
git checkout ../Biochemistry/compounds.* ../Biochemistry/Aliases/Unique_ModelSEED_Compound_*.txt

#Yeast Consensus Model (8.3.4)
../Scripts/Biochemistry/Add_New_Compounds.py Yeast_8.3.4_2019/yeastGEM_Compounds.tbl KEGG 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py Yeast_8.3.4_2019/yeastGEM_Reactions.tbl 'Published Model' f

#iJDZ836 (MetaCyc)
../Scripts/Biochemistry/Add_New_Compounds.py iJDZ836_23935467_2013/iJDZ836_Compounds.tbl MetaCyc 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iJDZ836_23935467_2013/iJDZ836_Reactions.tbl 'Published Model' f

#iAL1006 (RAVEN;InChI)
../Scripts/Biochemistry/Add_New_Compounds.py iAL1006_23555215_2013/iAL1006_Compounds.tbl InChI 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iAL1006_23555215_2013/iAL1006_Reactions.tbl 'Published Model' f

#iNX804 (RAVEN;KEGG)
../Scripts/Biochemistry/Add_New_Compounds.py iNX804_23172360_2013/iNX804_Compounds.tbl KEGG 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iNX804_23172360_2013/iNX804_Reactions.tbl 'Published Model' f

#iJB1325 (RAVEN;KEGG)
../Scripts/Biochemistry/Add_New_Compounds.py iJB1325_30275963_2018/iJB1325_ATCC1015_Compounds.tbl KEGG 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iJB1325_30275963_2018/iJB1325_ATCC1015_Reactions.tbl 'Published Model' f

#iWV1314 (RAVEN)
../Scripts/Biochemistry/Add_New_Compounds.py iWV1314_18500999_2008/iWV1314_Compounds.tbl iJB1325 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iWV1314_18500999_2008/iWV1314_Reactions.tbl 'Published Model' f

#iMA871 (RAVEN)
../Scripts/Biochemistry/Add_New_Compounds.py iMA871_18364712_2008/iMA871_Compounds.tbl iJB1325 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iMA871_18364712_2008/iMA871_Reactions.tbl 'Published Model' f

#iWV1213 (RAVEN)
../Scripts/Biochemistry/Add_New_Compounds.py iWV1213_26911256_2016/iWV1213_Compounds.tbl iJB1325 0 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iWV1213_26911256_2016/iWV1213_Reactions.tbl 'Published Model' f

#iYL619 (BiGG)
../Scripts/Biochemistry/Add_New_Compounds.py iYL619_PCP_23236514_2012/iYL619_PCP_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iYL619_PCP_23236514_2012/iYL619_PCP_Reactions.tbl 'Published Model' f

#iJL1454 (BiGG)
../Scripts/Biochemistry/Add_New_Compounds.py iJL1454_23624532_2013/iJL1454_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iJL1454_23624532_2013/iJL1454_Reactions.tbl 'Published Model' f

#iCY1106 (BiGG)
../Scripts/Biochemistry/Add_New_Compounds.py iCY1106_25582171_2015/iCY1106_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iCY1106_25582171_2015/iCY1106_Reactions.tbl 'Published Model' f

#iCT646 (BiGG)
../Scripts/Biochemistry/Add_New_Compounds.py iCT646_26915092_2016/iCT646_Compounds.tbl BiGG 1 'Published Model' f
../Scripts/Biochemistry/Add_New_Reactions.py iCT646_26915092_2016/iCT646_Reactions.tbl 'Published Model' f
