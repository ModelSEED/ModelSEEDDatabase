#!/bin/bash
#Copied from: https://gist.github.com/matsen/4263955
#Convert Main panel from SVG to TIFF
inkscape --without-gui --export-png=MSD_Thermodynamics_Rxn_dG.png --export-dpi 600 MSD_Thermodynamics_Rxn_dG.svg
convert -compress LZW -alpha remove MSD_Thermodynamics_Rxn_dG.png MSD_Thermodynamics_Rxn_dG.tiff
mogrify -alpha off MSD_Thermodynamics_Rxn_dG.tiff

#Convert inset panel from SVG to TIFF
inkscape --without-gui --export-png=MSD_Thermodynamics_Rxn_dGe.png --export-dpi 600 MSD_Thermodynamics_Rxn_dGe.svg
convert -compress LZW -alpha remove MSD_Thermodynamics_Rxn_dGe.png MSD_Thermodynamics_Rxn_dGe.tiff
mogrify -alpha off MSD_Thermodynamics_Rxn_dGe.tiff

#Build composite image
convert -resize 33% MSD_Thermodynamics_Rxn_dGe.tiff MSD_Thermodynamics_Rxn_dGe_Small.tiff
convert MSD_Thermodynamics_Rxn_dG.tiff MSD_Thermodynamics_Rxn_dGe_Small.tiff -gravity southeast -geometry +200+350 -composite MSD_Thermodynamics_Figure.tiff
