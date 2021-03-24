#!/usr/bin/env python
import sys

data=dict()
with open('../../Tables/Growth_Stats.tsv') as fh:
    for line in fh.readlines():
        line=line.strip('\r\n')
        array = line.split('\t')
        if(array[0] not in data):
            data[array[0]]=list()
        data[array[0]].append(array[2])
        data[array[0]].append(array[3])

categories = ['Compounds','Compounds with structure','Reactions','Mass-balanced reactions']
years = sorted(data.keys())

from bokeh.io import show, output_file, save, export_svgs
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.palettes import Colorblind, Set1

x = [ (cat, year) for cat in categories for year in years ]
counts = sum(zip(data[years[0]], data[years[1]], data[years[2]]), ()) # like an hstack

source = ColumnDataSource(data=dict(x=x, counts=counts))

p = figure(x_range=FactorRange(*x), plot_height=250, title=None, toolbar_location=None, tools="")

p.vbar(x='x', top='counts', width=0.9, source=source, line_color="white",

       # use the palette to colormap based on the the x[1:2] values
       fill_color=factor_cmap('x', palette=Colorblind[3], factors=years, start=1, end=2))

p.y_range.start = 0
p.x_range.range_padding = 0.1
p.xaxis.major_label_orientation = 1
p.xgrid.grid_line_color = None

file_name="MSD_Growth"
output_file(file_name+".html")
save(p)

p.toolbar.logo = None
p.toolbar_location = None

p.background_fill_color = None
p.border_fill_color = None

p.output_backend="svg"
export_svgs(p,filename=file_name+".svg")
