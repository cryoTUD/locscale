#!/bin/bash
echo "Adding proper path to your bashrc file"
locscale=$(pwd)
echo $locscale
echo 'export PATH=$PATH:'$locscale >> ~/.bashrc
echo 'export PYTHONPATH=${PYTHONPATH}:'$locscale >> ~/.bashrc
scripts=$locscale/scripts
utils=$scripts/utils
pseudomodel=$scripts/get_pseudomodel
echo 'export PATH=${PATH}:'$scripts >> ~/.bashrc
echo 'export PATH=${PATH}:'$utils >> ~/.bashrc
echo 'export PATH=${PATH}:'$pseudomodel >> ~/.bashrc
echo 'export PYTHONPATH=${PYTHONPATH}:'$scripts >> ~/.bashrc
echo 'export PYTHONPATH=${PYTHONPATH}:'$utils >> ~/.bashrc
echo 'export PYTHONPATH=${PYTHONPATH}:'$pseudomodel >> ~/.bashrc

