#!/bin/bash
#Copied from: https://gist.github.com/matsen/4263955
#Convert Main panel from SVG to TIFF
inkscape --without-gui --export-png=MSD_Growth.png --export-dpi 600 MSD_Growth.svg
convert -compress LZW -alpha remove MSD_Growth.png MSD_Growth.tiff
mogrify -alpha off MSD_Growth.tiff
