import argparse
sample_run_locscale = "python /path/to/locscale/main.py run_locscale --emmap_path path/to/emmap.mrc -res 3.4 -o locscale.mrc --verbose"
sample_run_emmernet = "python /path/to/locscale/main.py run_emmernet --emmap_path path/to/emmap.mrc --verbose"
description = ["*** Optimisation of contrast in cryo-EM density maps using local density scaling ***\n",\
    "Command line arguments: \n",\
        "LocScale: \n",\
        "{}\n".format(sample_run_locscale),\
        "EMmerNet: \n",\
        "{}".format(sample_run_emmernet)]

main_parser = argparse.ArgumentParser(prog="locscale",
description="".join(description)) 

## Add subparsers
sub_parser = main_parser.add_subparsers(dest='command')
locscale_parser = sub_parser.add_parser('run_locscale', help='Run LocScale')
emmernet_parser = sub_parser.add_parser('run_emmernet', help='Run EMMERNET')
test_parser = sub_parser.add_parser('test', help='Run tests')

# **************************************************************************************
# ************************ Command line arguments LocScale *****************************
# **************************************************************************************

## Input either unsharpened EM map or two halfmaps
locscale_emmap_input = locscale_parser.add_mutually_exclusive_group(required=True)
locscale_emmap_input.add_argument('-em', '--emmap_path',  help='Path to unsharpened EM map')
locscale_emmap_input.add_argument('-hm', '--halfmap_paths', nargs=2, help='Paths to first and second halfmaps')

## Input model map file (mrc file) or atomic model (pdb file)
locscale_parser.add_argument('-mm', '--model_map', help='Path to model map file')
locscale_parser.add_argument('-mc', '--model_coordinates', help='Path to PDB file', default=None)

## Input mask
locscale_parser.add_argument('-ma', '--mask', help='Input filename mask')

## Output arguments
locscale_parser.add_argument('-o', '--outfile', help='Output filename', default="locscale_output.mrc")
locscale_parser.add_argument('-v', '--verbose', action='store_true',help='Verbose output')
locscale_parser.add_argument('--report_filename', type=str, help='Filename for storing PDF output and statistics', default="locscale_report")
locscale_parser.add_argument('-op', '--output_processing_files', type=str, help='Path to store processing files', default=None)

## LocScale main function parameters
locscale_parser.add_argument('-wn', '--window_size', type=int, help='window size in pixels', default=None)
locscale_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,help='MPI version')
locscale_parser.add_argument('-np', '--number_processes', help='Number of processes to use', type=int, default=1)

## Refinement parameters
locscale_parser.add_argument('-ref_it', '--refmac_iterations', help='For atomic model refinement: number of refmac iterations', default=10, type=int)
locscale_parser.add_argument('-res', '--ref_resolution', type=float, help='Resolution target for Refmac refinement')
locscale_parser.add_argument('-p', '--apix', type=float, help='pixel size in Angstrom')
locscale_parser.add_argument('--add_blur', type=int, help='Globally sharpen the target map for REFMAC refinement', default=20)
locscale_parser.add_argument('--refmac5_path', type=str, help='Path to refmac5 executable', default=None)
locscale_parser.add_argument('--cref_pickle', type=str, help='Path for Cref filter for the target map of bfactor refinement', default=None)


## Model map parameters
locscale_parser.add_argument('-mres', '--model_resolution', type=float, help='Resolution limit for Model Map generation')
locscale_parser.add_argument('-sym', '--symmetry', default='C1', type=str, help='Impose symmetry condition for output')

## FDR parameters
locscale_parser.add_argument('-fdr_w', '--fdr_window_size', type=int, help='window size in pixels for FDR thresholding', default=None)
locscale_parser.add_argument('--averaging_filter_size', '--averaging_filter_size', type=int, help='window size in pixels for FDR thresholding', default=3)
locscale_parser.add_argument('-fdr_f', '--fdr_filter', type=float, help='Pre-filter for FDR thresholding', default=None)
locscale_parser.add_argument('-th', '--mask_threshold', type=float, help='FDR threshold', default=0.5)

## Integrated pseudo-atomic model method parameters
locscale_parser.add_argument('--complete_model', help='Add pseudo-atoms to areas of the map which are not modelled', action='store_true')
locscale_parser.add_argument('-avg_w', '--averaging_window', type=int, help='Window size for filtering the fdr difference map for integrated pseudo-model', default=3)

## Pseudo-atomic model method parameters
locscale_parser.add_argument('-pm', '--pseudomodel_method', help='For pseudo-atomic model: method', default='gradient')
locscale_parser.add_argument('-pm_it', '--total_iterations', type=int, help='For pseudo-atomic model: total iterations', default=50)
locscale_parser.add_argument('-dst', '--distance', type=float, help='For pseudo-atomic model: typical distance between atoms', default=1.2)
locscale_parser.add_argument('-mw', '--molecular_weight', help='Input molecular weight (in kDa)', default=None, type=float)
locscale_parser.add_argument('--build_ca_only', help='For gradient pseudomodel building: use only Ca atoms with interatomic distance 3.8', action='store_true',default=False)
locscale_parser.add_argument('-s', '--smooth_factor', type=float, help='Smooth factor for merging profiles', default=0.3)
locscale_parser.add_argument('--boost_secondary_structure', type=float, help='Amplify signal corresponding to secondary structures', default=1.5)
locscale_parser.add_argument('--no_reference', action='store_true', default=False,help='Run locscale without using any reference information')
locscale_parser.add_argument('--set_local_bfactor', type=float, default=20,help='For reference-less sharpening. Use this value to set the local b-factor of the maps')

## non-default arguments
locscale_parser.add_argument('--dev_mode', action='store_true', default=False,help='If true, this will force locscale to use the theoretical profile even if model map present and will not check for user input consistency')
locscale_parser.add_argument('--skip_refine', help='Ignore REFMAC refinement', action='store_true')

# **************************************************************************************
# ************************ Command line arguments EMMERNET *******************************
# **************************************************************************************

## Input either unsharpened EM map or two halfmaps
emmernet_emmap_input = emmernet_parser.add_mutually_exclusive_group(required=True)
emmernet_emmap_input.add_argument('-em', '--emmap_path',  help='Path to unsharpened EM map')
emmernet_emmap_input.add_argument('-hm', '--halfmap_paths', nargs=2, help='Paths to first and second halfmaps')

## Output arguments
emmernet_parser.add_argument('-o', '--outfile', help='Output filename', default="emmernet_output.mrc")
emmernet_parser.add_argument('-op', '--output_processing_files', type=str, help='Path to store processing files', default=None)
emmernet_parser.add_argument('-v', '--verbose', action='store_true',help='Verbose output')

## Emmernet main function parameters
emmernet_parser.add_argument('-trained_model','--trained_model', help='Type of emmernet model to use', \
                            choices=['hybrid','model_based', 'model_free', 'ensemble','model_based_no_freqaug'], default='hybrid')
emmernet_parser.add_argument('-s', '--stride', help='Stride for EMMERNET', default=16, type=int)
emmernet_parser.add_argument('-bs', '--batch_size', type=int, help='Batch size for EMMERNET', default=8)
emmernet_parser.add_argument("-gpus", "--gpu_ids", nargs='+', help="numbers of the selected GPUs, format: '1 2 3 ... 5'", required=False)
emmernet_parser.add_argument('-target', '--target_map_path', type=str, help='Path to the target map for phase correlations', default=None)
emmernet_parser.add_argument('-download', '--download', help='Download the model weights', action='store_true', default=False)

############################################################################################
# ************************ Command line arguments TESTS ********************************** #
############################################################################################

