#!/bin/bash
MODEL=${1}
NAME=${2}
MAP=${3}
RES=${4}
CCPEM_PATH=${5}
CCP4_PATH=${6}
H=${7}
W=${8}
L=${9}
NUM_ITER=${10}
#Step 1: 
REFMAC_SCRIPT=${CCP4_PATH}"/bin/refmac5"
echo ${CCP4_PATH}
echo ${REFMAC_SCRIPT}
$REFMAC_SCRIPT mapin $MAP hklout starting_map.mtz  << eof
mode sfcalc
source EM MB
reso ${RES}
sfcalc blur
end
eof
#Step 2:
${CCP4_PATH}/bin/pdbset xyzin $MODEL xyzout pdbset.pdb  << eof
CELL  ${H} ${W} ${L} 90.0 90.0 90.0
end
eof
REFMAC_OUTPUT="${NAME}_refmac_refined.pdb"
#Step 3:
${CCP4_PATH}/bin/refmac5 xyzin pdbset.pdb hklin starting_map.mtz xyzout $REFMAC_OUTPUT  << eof
labin FP=Fout0 PHIB=Pout0
make hydr no
solvent no
source EM MB
ncycle ${NUM_ITER}
weight auto
ncsr local
reso ${RES}
ridge dist sigma 0.01
ridge dist dmax 4.2
BFACtor SET 40
end
eof

