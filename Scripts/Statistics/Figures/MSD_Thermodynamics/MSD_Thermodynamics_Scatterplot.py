#!/usr/bin/env python
import sys
from BiochemPy import Reactions
reactions_helper=Reactions()
reactions_dict=reactions_helper.loadReactions()

thermo_root="../../../Thermodynamics"
output_file_root="MSD_Thermodynamics"
type_revs_rxns=dict()
with open(thermo_root+"/Estimated_Reaction_Reversibility_Report_EQ.txt") as fh:
    for line in fh.readlines():
        line=line.strip()
        (rxn,rule,gf_rev,eq_rev) = line.split('\t')
        type_rev = gf_rev+eq_rev

        if("OK" not in reactions_dict[rxn]['status']):
            continue

        if(reactions_dict[rxn]['is_obsolete']==1):
            continue

        if(reactions_dict[rxn]['is_transport']==1):
            continue

        if('GFC' not in reactions_dict[rxn]['notes'] or 'EQU' not in reactions_dict[rxn]['notes']):
            continue

        if(gf_rev==eq_rev):
            type_rev="no-change"
        else:
            type_rev="irreversible-to-irreversible"
            if(eq_rev=="="):
                type_rev="irreverible-to-reversible"
            if(gf_rev=="="):
                type_rev="reversible-to-irreversible"

        if(type_rev not in type_revs_rxns):
            type_revs_rxns[type_rev]=list()
        type_revs_rxns[type_rev].append(rxn)
fh.close()

data_categories=dict()
rxn_type_cpds=dict()
(max_dg,min_dg)=(-10000,10000)
(max_dge,min_dge)=(-10000,10000)

with open(thermo_root+"/Reactions_GroupFormation_eQuilibrator_Comparison.txt") as fh:
    header=1
    for line in fh.readlines():
        if(header==1):
            header=0
            continue

        line=line.strip()
        (rxn,gf_en,eq_en)=line.split('\t')

        if(gf_en=="nan" or eq_en=="nan"):
            continue

        rev_type=None
        for key in type_revs_rxns:
            if(rxn in type_revs_rxns[key]):
                rev_type = key

        if(rev_type is None):
            continue

        if(rev_type not in rxn_type_cpds):
            rxn_type_cpds[rev_type]=dict()
        rxn_type_cpds[rev_type][rxn]=reactions_dict[rxn]['compound_ids'].split(';')

        (gf_dg,gf_dge) = gf_en.split('|')
        (eq_dg,eq_dge) = eq_en.split('|')

        if(rev_type  not in data_categories):
            data_categories[rev_type]={'dg':{'x':[],'y':[]},'dge':{'x':[],'y':[]}}

        data_categories[rev_type]['dg']['x'].append(float(eq_dg))
        data_categories[rev_type]['dg']['y'].append(float(gf_dg))
        data_categories[rev_type]['dge']['x'].append(float(eq_dge))
        data_categories[rev_type]['dge']['y'].append(float(gf_dge))

        if(float(eq_dg) > max_dg):
            max_dg=float(eq_dg)
        if(float(eq_dge) > max_dge):
            max_dge=float(eq_dge)
        if(float(gf_dg) > max_dg):
            max_dg=float(gf_dg)
        if(float(gf_dge) > max_dge):
            max_dge=float(gf_dge)

        if(float(eq_dg) < min_dg):
            min_dg=float(eq_dg)
        if(float(eq_dge) < min_dge):
            min_dge=float(eq_dge)
        if(float(gf_dg) < min_dg):
            min_dg=float(gf_dg)
        if(float(gf_dge) < min_dge):
            min_dge=float(gf_dge)
fh.close()

print("Rxn dG Range: "+str(min_dg)+" - "+str(max_dg))
print("Rxn dGe Range: "+str(min_dge)+" - "+str(max_dge))

from bokeh.plotting import figure, output_file, show, save
from bokeh.io import export_svgs, export_png
from bokeh.palettes import Colorblind

color_dict={'no-change':Colorblind[4][0],
            'irreversible-to-irreversible':Colorblind[4][1],
            'irreverible-to-reversible':Colorblind[4][2],
            'reversible-to-irreversible':Colorblind[4][3]}

#Rxn dG
dg_plot = figure(plot_width=500, plot_height=500,x_range=(-500,500),y_range=(-500,500))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dg']['x']), reverse=True)):
    dg_plot.circle( data_categories[category]['dg']['x'], data_categories[category]['dg']['y'],
                    color=color_dict[category],legend=category+" ("+str(len(data_categories[category]['dg']['x']))+")")

dg_plot.xaxis.axis_label = 'eQuilibrator dG'
dg_plot.yaxis.axis_label = 'Group Formation dG'

dg_plot.legend.location = "top_left"
dg_plot.legend.click_policy="hide"

file_name=output_file_root+"_Rxn_dG"
output_file(file_name+".html")
save(dg_plot)

dg_plot.toolbar.logo = None
dg_plot.toolbar_location = None

dg_plot.background_fill_color = None
dg_plot.border_fill_color = None

export_png(dg_plot,filename=file_name+".png")
dg_plot.output_backend="svg"
export_svgs(dg_plot,filename=file_name+".svg")

#Rxn dGe
dge_plot = figure(plot_width=500, plot_height=500,x_range=(0,120),y_range=(0,120))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dge']['x']), reverse=True)):
    dge_plot.circle( data_categories[category]['dge']['x'], data_categories[category]['dge']['y'],
                     color=color_dict[category],legend=category+" ("+str(len(data_categories[category]['dge']['x']))+")")

dge_plot.xaxis.axis_label = 'eQuilibrator dGe'
dge_plot.yaxis.axis_label = 'Group Formation dGe'

dge_plot.legend.location = "top_right"
dge_plot.legend.click_policy="hide"

file_name=output_file_root+"_Rxn_dGe"
output_file(file_name+".html")
save(dge_plot)

dge_plot.toolbar.logo = None
dge_plot.toolbar_location = None

dge_plot.background_fill_color = None
dge_plot.border_fill_color = None

export_png(dge_plot,filename=file_name+".png")
dge_plot.output_backend="svg"
export_svgs(dge_plot,filename=file_name+".svg")

from BiochemPy import Compounds
compounds_helper = Compounds()
compounds_dict=compounds_helper.loadCompounds()

compounds_energies=dict()
with open(thermo_root+"/Compounds_GroupFormation_eQuilibrator_Comparison.txt") as fh:
    header=1
    for line in fh.readlines():
        if(header==1):
            header=0
            continue

        line=line.strip()
        (cpd,gf_en,eq_en)=line.split('\t')

        if(gf_en=="nan" or eq_en=="nan"):
            continue

        (gf_dg,gf_dge) = gf_en.split('|')
        (eq_dg,eq_dge) = eq_en.split('|')
        compounds_energies[cpd]={'gf':[float(gf_dg),float(gf_dge)],'eq':[float(eq_dg),float(eq_dge)]}

sys.exit()
categories=data_categories.keys()
data_categories=dict()
touched_cpds=dict()
(max_dg,min_dg)=(-10000,10000)
(max_dge,min_dge)=(-10000,10000)

for rev_type in categories:
    if(rev_type not in data_categories):
        data_categories[rev_type]={'dg':{'x':[],'y':[]},'dge':{'x':[],'y':[]}}

    for rxn in rxn_type_cpds[rev_type]:
        for cpd in rxn_type_cpds[rev_type][rxn]:
            if(cpd not in compounds_energies or cpd in touched_cpds):
                continue
            touched_cpds[cpd]=1

            if(compounds_energies[cpd]['eq'][0] > max_dg):
                max_dg=compounds_energies[cpd]['eq'][0]
            if(compounds_energies[cpd]['eq'][1] > max_dge):
                max_dge=compounds_energies[cpd]['eq'][1]
            if(compounds_energies[cpd]['gf'][0] > max_dg):
                max_dg=compounds_energies[cpd]['gf'][0]
            if(compounds_energies[cpd]['gf'][1] > max_dge):
                max_dge=compounds_energies[cpd]['gf'][1]

            if(compounds_energies[cpd]['eq'][0] < min_dg):
                min_dg=compounds_energies[cpd]['eq'][0]
            if(compounds_energies[cpd]['eq'][1] < min_dge):
                min_dge=compounds_energies[cpd]['eq'][1]
            if(compounds_energies[cpd]['gf'][0] < min_dg):
                min_dg=compounds_energies[cpd]['gf'][0]
            if(compounds_energies[cpd]['gf'][1] < min_dge):
                min_dge=compounds_energies[cpd]['gf'][1]

            data_categories[rev_type]['dg']['x'].append(compounds_energies[cpd]['eq'][0])
            data_categories[rev_type]['dg']['y'].append(compounds_energies[cpd]['gf'][0])
            data_categories[rev_type]['dge']['x'].append(compounds_energies[cpd]['eq'][1])
            data_categories[rev_type]['dge']['y'].append(compounds_energies[cpd]['gf'][1])

print("Cpd dG Range: "+str(min_dg)+" - "+str(max_dg))
print("Cpd dGe Range: "+str(min_dge)+" - "+str(max_dge))

#Cpd dG
dg_plot = figure(plot_width=500, plot_height=500,x_range=(-3200,3200),y_range=(-3200,3200))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dg']['x']), reverse=True)):
    dg_plot.circle( data_categories[category]['dg']['x'], data_categories[category]['dg']['y'],
                    color=color_dict[category],legend=category+" ("+str(len(data_categories[category]['dg']['x']))+")")

dg_plot.xaxis.axis_label = 'eQuilibrator dG'
dg_plot.yaxis.axis_label = 'Group Formation dG'

dg_plot.legend.location = "top_left"
dg_plot.legend.click_policy="hide"

file_name=output_file_root+"_Cpd_dG"
output_file(file_name+".html")
save(dg_plot)

dg_plot.toolbar.logo = None
dg_plot.toolbar_location = None

dg_plot.background_fill_color = None
dg_plot.border_fill_color = None

export_png(dg_plot,filename=file_name+".png")
dg_plot.output_backend="svg"
export_svgs(dg_plot,filename=file_name+".svg")

#Cpd dGe
dge_plot = figure(plot_width=500, plot_height=500,x_range=(0,100000),y_range=(0,100000))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dge']['x']), reverse=True)):
    dge_plot.circle( data_categories[category]['dge']['x'], data_categories[category]['dge']['y'], 
                     color=color_dict[category],legend=category+" ("+str(len(data_categories[category]['dge']['x']))+")")

dge_plot.xaxis.axis_label = 'eQuilibrator dGe'
dge_plot.yaxis.axis_label = 'Group Formation dGe'

dge_plot.legend.location = "top_right"
dge_plot.legend.click_policy="hide"

file_name=output_file_root+"_Cpd_dGe"
output_file(file_name+".html")
save(dge_plot)

dge_plot.toolbar.logo = None
dge_plot.toolbar_location = None

dge_plot.background_fill_color = None
dge_plot.border_fill_color = None

export_png(dge_plot,filename=file_name+".png")
dge_plot.output_backend="svg"
export_svgs(dge_plot,filename=file_name+".svg")
