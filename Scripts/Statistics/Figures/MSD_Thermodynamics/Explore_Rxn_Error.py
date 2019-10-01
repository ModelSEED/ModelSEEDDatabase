#!/usr/bin/env python
thermo_root="../../../Thermodynamics"
with open(thermo_root+"/Reactions_GroupFormation_eQuilibrator_Comparison.txt") as fh:
    header=1
    for line in fh.readlines():
        if(header==1):
            header=0
            continue

        line=line.strip()
        (cpd,gf_en,eq_en)=line.split('\t')

        if(gf_en != "nan" and "10000000.0" not in gf_en):
            (gf_dg,gf_dge) = gf_en.split('|')
            if(float(gf_dg)!=0):
                print(float(gf_dge)/abs(float(gf_dg)),gf_dg,gf_dge,"GF")

        if(eq_en != "nan"):
            (eq_dg,eq_dge) = eq_en.split('|')
            if(float(eq_dg)!=0):
                print(float(eq_dge)/abs(float(eq_dg)),eq_dg,eq_dge,"EQ")

        
        if(gf_en != "nan" and "10000000.0" not in gf_en and eq_en != "nan"):
            print(abs(float(gf_dg)-float(eq_dg)),"DIFF")
