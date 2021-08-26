#!/bin/bash
map1=$1
map2=$2
output=$3
apix=$4
echo "Script to add two half maps. Usage: add_half_maps.sh halfmap1.mrc halfmap2.mrc outputname.mrc <apix-in-float>"
e2proc3d.py $map1 $output --addfile=$map2 --apix=$apix

