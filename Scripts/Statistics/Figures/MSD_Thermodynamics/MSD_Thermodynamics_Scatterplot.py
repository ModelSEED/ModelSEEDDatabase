#!/usr/bin/env python
import sys
from BiochemPy import Reactions
reactions_helper=Reactions()
reactions_dict=reactions_helper.loadReactions()

thermo_root="../../../Thermodynamics"
output_file_root="MSD_Thermodynamics"
type_revs_rxns=dict()
accepted_rxns=dict()
with open(thermo_root+"/Estimated_Reaction_Reversibility_Report_EQ.txt") as fh:
    for line in fh.readlines():
        line=line.strip()
        (rxn,rule,gc_rev,eq_rev) = line.split('\t')
        type_rev = gc_rev+eq_rev

        if(rxn not in reactions_dict):
            continue

        if(reactions_dict[rxn]['status'] == "EMPTY"):
            continue

        if("OK" not in reactions_dict[rxn]['status']):
            continue

        if(reactions_dict[rxn]['is_obsolete']==1):
            continue

        if(reactions_dict[rxn]['is_transport']==1):
            continue

        if('GCC' not in reactions_dict[rxn]['notes'] or 'EQU' not in reactions_dict[rxn]['notes']):
            continue

        accepted_rxns[rxn]=1

        if(gc_rev==eq_rev):
            type_rev="no-change"
        else:
            type_rev="irreversible-to-irreversible"
            if(eq_rev=="="):
                type_rev="irreverible-to-reversible"
            if(gc_rev=="="):
                type_rev="reversible-to-irreversible"

        if(type_rev not in type_revs_rxns):
            type_revs_rxns[type_rev]=list()
        type_revs_rxns[type_rev].append(rxn)
fh.close()

data_categories=dict()
rxn_type_cpds=dict()

min_max_dict={'gc_dg':{'max':-sys.maxsize-1,'min':sys.maxsize},
              'eq_dg':{'max':-sys.maxsize-1,'min':sys.maxsize},
              'gc_dge':{'max':-sys.maxsize-1,'min':sys.maxsize},
              'eq_dge':{'max':-sys.maxsize-1,'min':sys.maxsize},
              'gc_dge_pct':{'max':-sys.maxsize-1,'min':sys.maxsize},
              'eq_dge_pct':{'max':-sys.maxsize-1,'min':sys.maxsize}}

total_count=0
included_count=0

diff_dict={'small':[],'medium':[],'large':[]}
high_pct={'gc':[],'eq':[]}
with open(thermo_root+"/Reactions_GroupContribution_eQuilibrator_Comparison.txt") as fh:
    header=1
    for line in fh.readlines():
        if(header==1):
            header=0
            continue

        line=line.strip()
        (rxn,gc_en,eq_en)=line.split('\t')

        if(rxn not in accepted_rxns):
            continue

        if(gc_en=="nan" or eq_en=="nan"):
            del accepted_rxns[rxn]
            continue

        rev_type=None
        for key in type_revs_rxns:
            if(rxn in type_revs_rxns[key]):
                rev_type = key

        if(rev_type is None):
            print("Warning: "+rxn+" missing")
            continue

        if(rev_type not in rxn_type_cpds):
            rxn_type_cpds[rev_type]=dict()
        rxn_type_cpds[rev_type][rxn]=reactions_dict[rxn]['compound_ids'].split(';')

        (gc_dg,gc_dge) = gc_en.split('|')
        (eq_dg,eq_dge) = eq_en.split('|')        

        if(rev_type  not in data_categories):
            data_categories[rev_type]={'dg':{'x':[],'y':[]},'dge':{'x':[],'y':[]}}

        data_categories[rev_type]['dg']['x'].append(float(eq_dg))
        data_categories[rev_type]['dg']['y'].append(float(gc_dg))

        #Convert to Percentage
        gc_dge_pct = 0.0
        eq_dge_pct = 0.0

        #If the value is zero, I'm assuming a small number as denominator
        if(abs( float(gc_dg) ) > 0):
            gc_dge_pct = ( float(gc_dge)/ abs(float(gc_dg)) ) * 100.0
        else:
            gc_dge_pct = ( float(gc_dge) / 0.001 ) * 100.0

        if(abs( float(eq_dg) ) > 0):
            eq_dge_pct = ( float(eq_dge)/ abs(float(eq_dg)) ) * 100.0
        else:
            eq_dge_pct = ( float(eq_dge) / 0.001 ) * 100.0

        total_count+=1
        if(eq_dge_pct < 10000 and gc_dge_pct < 10000):
            included_count+=1

        data_categories[rev_type]['dge']['x'].append(float(eq_dge_pct))
        data_categories[rev_type]['dge']['y'].append(float(gc_dge_pct))

        if(float(gc_dg) > min_max_dict['gc_dg']['max']):
            min_max_dict['gc_dg']['max']=float(gc_dg)
        if(float(eq_dg) > min_max_dict['eq_dg']['max']):
            min_max_dict['eq_dg']['max']=float(eq_dg)
        if(float(gc_dge) > min_max_dict['gc_dge']['max']):
            min_max_dict['gc_dge']['max']=float(gc_dge)
        if(float(eq_dge) > min_max_dict['eq_dge']['max']):
            min_max_dict['eq_dge']['max']=float(eq_dge)
        if(float(gc_dge_pct) > min_max_dict['gc_dge_pct']['max']):
            min_max_dict['gc_dge_pct']['max']=float(gc_dge_pct)
        if(float(eq_dge_pct) > min_max_dict['eq_dge_pct']['max']):
            min_max_dict['eq_dge_pct']['max']=float(eq_dge_pct)

        if(float(gc_dg) < min_max_dict['gc_dg']['min']):
            min_max_dict['gc_dg']['min']=float(gc_dg)
        if(float(eq_dg) < min_max_dict['eq_dg']['min']):
            min_max_dict['eq_dg']['min']=float(eq_dg)
        if(float(gc_dge) < min_max_dict['gc_dge']['min']):
            min_max_dict['gc_dge']['min']=float(gc_dge)
        if(float(eq_dge) < min_max_dict['eq_dge']['min']):
            min_max_dict['eq_dge']['min']=float(eq_dge)
        if(float(gc_dge_pct) < min_max_dict['gc_dge_pct']['min']):
            min_max_dict['gc_dge_pct']['min']=float(gc_dge_pct)
        if(float(eq_dge_pct) < min_max_dict['eq_dge_pct']['min']):
            min_max_dict['eq_dge_pct']['min']=float(eq_dge_pct)

        gc_eq_diff = abs( float(gc_dg) - float(eq_dg) )
        if(gc_eq_diff < 5):
            diff_dict['small'].append(eq_dge_pct)
        elif(gc_eq_diff < 15):
            diff_dict['medium'].append(eq_dge_pct)
        elif(gc_eq_diff >= 15):
            diff_dict['large'].append(eq_dge_pct)

        if(gc_dge_pct > 100.0):
            high_pct['gc'].append(gc_dge_pct)
        if(eq_dge_pct > 100.0):
            high_pct['eq'].append(eq_dge_pct)

fh.close()

print("Fraction in dGe panel: "+str(float(included_count)/float(total_count)))
print(str(len(accepted_rxns)))
for size_diff in diff_dict:
    print(size_diff,str(len(diff_dict[size_diff])))

large_uncertainty_count=0
for pct in diff_dict['large']:
    if(pct > 100):
        large_uncertainty_count+=1
print(large_uncertainty_count)

print(str(len(high_pct['gc'])),str(len(high_pct['eq'])))
print("Rxn Ranges: ",min_max_dict)

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dg']['x']), reverse=True)):
    print(category,str(len(data_categories[category]['dg']['x'])))

from bokeh.plotting import figure, output_file, show, save
from bokeh.io import export_svgs
from bokeh.palettes import Colorblind

color_dict={'no-change':Colorblind[4][0],
            'irreversible-to-irreversible':Colorblind[4][1],
            'irreverible-to-reversible':Colorblind[4][2],
            'reversible-to-irreversible':Colorblind[4][3]}

#Rxn dG
dg_plot = figure(plot_width=500, plot_height=500,x_range=(-500,500),y_range=(-500,500))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dg']['x']), reverse=True)):
    dg_plot.circle( data_categories[category]['dg']['x'], data_categories[category]['dg']['y'],
                    color=color_dict[category],legend_label=category+" ("+str(len(data_categories[category]['dg']['x']))+")")

# UniCode names taken from here: https://unicode.org/Public/UNIDATA/NamesList.txt
deltag_string = '\N{greek capital letter delta}\N{latin subscript small letter r}G\N{prime}\N{superscript zero}'
dg_plot.xaxis.axis_label = deltag_string+' (eQuilibrator)'
dg_plot.yaxis.axis_label = deltag_string+' (Group Contribution)'

dg_plot.legend.location = "top_left"
dg_plot.legend.click_policy="hide"

file_name=output_file_root+"_Rxn_dG"
output_file(file_name+".html")
save(dg_plot)

dg_plot.toolbar.logo = None
dg_plot.toolbar_location = None

dg_plot.background_fill_color = None
dg_plot.border_fill_color = None

dg_plot.output_backend="svg"
export_svgs(dg_plot,filename=file_name+".svg")

#Rxn dGe
dge_plot = figure(plot_width=500, plot_height=500,x_range=(-100,10000),y_range=(-100,10000))

for i, category in enumerate(sorted(data_categories, key=lambda k: len(data_categories[k]['dge']['x']), reverse=True)):
    dge_plot.circle( data_categories[category]['dge']['x'], data_categories[category]['dge']['y'],
                     color=color_dict[category] ) #,legend_label=category+" ("+str(len(data_categories[category]['dge']['x']))+")")

#dge_plot.xaxis.axis_label = deltag_string+' % uncertainty (eQuilibrator)'
#dge_plot.yaxis.axis_label = deltag_string+' % uncertainty (Group Contribution)'
dge_plot.xaxis.axis_label = '% Uncertainty (eQuilibrator)'
dge_plot.yaxis.axis_label = '% Uncertainty (Group Contribution)'
dge_plot.xaxis.axis_label_text_font_size = "14pt"
dge_plot.yaxis.axis_label_text_font_size = "14pt"

#dge_plot.legend.location = "top_right"
#dge_plot.legend.click_policy="hide"

file_name=output_file_root+"_Rxn_dGe"
output_file(file_name+".html")
save(dge_plot)

dge_plot.toolbar.logo = None
dge_plot.toolbar_location = None

dge_plot.background_fill_color = None
dge_plot.border_fill_color = None

dge_plot.output_backend="svg"
export_svgs(dge_plot,filename=file_name+".svg")

###############################################################
## Compounds Figure not being included in paper at this time
###############################################################
sys.exit()

from BiochemPy import Compounds
compounds_helper = Compounds()
compounds_dict=compounds_helper.loadCompounds()

compounds_energies=dict()
with open(thermo_root+"/Compounds_GroupContribution_eQuilibrator_Comparison.txt") as fh:
    header=1
    for line in fh.readlines():
        if(header==1):
            header=0
            continue

        line=line.strip()
        (cpd,gc_en,eq_en)=line.split('\t')

        if(compounds_dict[cpd]['is_obsolete']==1):
            continue

        if(gc_en=="nan" or eq_en=="nan"):
            continue

        (gc_dg,gc_dge) = gc_en.split('|')
        (eq_dg,eq_dge) = eq_en.split('|')
        compounds_energies[cpd]={'gf':[float(gc_dg),float(gc_dge)],'eq':[float(eq_dg),float(eq_dge)]}

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
                    color=color_dict[category],legend_label=category+" ("+str(len(data_categories[category]['dg']['x']))+")")

dg_plot.xaxis.axis_label = 'eQuilibrator dG'
dg_plot.yaxis.axis_label = 'Group Contribution dG'

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
                     color=color_dict[category],legend_label=category+" ("+str(len(data_categories[category]['dge']['x']))+")")

dge_plot.xaxis.axis_label = 'eQuilibrator dGe'
dge_plot.yaxis.axis_label = 'Group Contribution dGe'

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
