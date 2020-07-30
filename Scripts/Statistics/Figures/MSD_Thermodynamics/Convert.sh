#!/bin/bash
convert MSD_Thermodynamics_Rxn_dG.png -resize 2000x2000 -density 300 MSD_Thermodynamics_Rxn_dG.tiff
convert MSD_Thermodynamics_Rxn_dGe.png -resize 2000x2000 -density 300 MSD_Thermodynamics_Rxn_dGe.tiff
convert -resize 33% MSD_Thermodynamics_Rxn_dGe.tiff MSD_Thermodynamics_Rxn_dGe_Small.tiff
convert MSD_Thermodynamics_Rxn_dG.tiff MSD_Thermodynamics_Rxn_dGe_Small.tiff -gravity southeast -geometry +100+250 -composite output.tiff
