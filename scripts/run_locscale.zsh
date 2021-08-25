#!/bin/zsh
#Step1: Map Mask <REFMAC_NAME> 
echo "Locscale step 1: Map Mask"
CCP4_PATH=$8
CCPEM_PATH=$7
WSIZE=$6
PIX=$5
##Input_step1
REFMAP_NAME=$3
REFMAP_XYZ=$4
${CCP4_PATH}/bin/mapmask mapin1 $REFMAP_NAME mapout $REFMAP_XYZ  << eof
axis X Y Z
end
eof
EMMAP=$1
EMMAP_XYZ=$2
${CCP4_PATH}/bin/mapmask mapin1 $EMMAP mapout $EMMAP_XYZ  << eof
axis X Y Z
end
eof
##Step2: Locscale?
echo "Locscale step 2: Locscale"
${CCPEM_PATH}/bin/ccpem-mpirun -np 4 ${CCPEM_PATH}/bin/ccpem-python ${CCPEM_PATH}/lib/py2/loc_scale/np_locscale.py --em_map $EMMAP_XYZ --model_map $REFMAP_XYZ --apix $PIX --window_size $WSIZE --outfile loc_scale.mrc --verbose True --mpi 
