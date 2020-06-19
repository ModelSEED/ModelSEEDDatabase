To run the bokeh scripts and export to png and svg, you need these packages
```
pip install bokeh
conda install -c bokeh selenium
conda install -c conda-forge phantomjs
conda install -c conda-forge firefox geckodriver
```

(the last installation package was needed for selenium apparently)

Running MSD_Growth/MSD_Growth_BarChart.py will re-generate Figure 1
Running MSD_Thermodynamics/MSD_Thermodynamics_Scatterplot.py

The latter also prints out numbers we use for Figure 2 and the "Thermodynamics" section