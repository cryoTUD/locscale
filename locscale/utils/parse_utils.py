import argparse
sample_run_locscale = "python /path/to/locscale/main.py run_locscale --emmap_path path/to/emmap.mrc -res 3.4 -o locscale.mrc --verbose"
sample_run_emmernet = "python /path/to/locscale/main.py run_emmernet --emmap_path path/to/emmap.mrc --verbose"
description = ["*** Optimisation of contrast in cryo-EM density maps using local density scaling ***\n",\
    "Command line arguments: \n",\
        "LocScale: \n",\
        "{}\n".format(sample_run_locscale),\
        "EMmerNet: \n",\
        "{}".format(sample_run_emmernet)]

# **************************************************************************************
# ************************ Command line arguments LocScale *****************************
# **************************************************************************************
def add_common_arguments(parser):
    ## Input either unsharpened EM map or two halfmaps
    input_emmap_group = parser.add_argument_group('Map input arguments')
    locscale_emmap_input = input_emmap_group.add_mutually_exclusive_group(required=False)
    locscale_emmap_input.add_argument(
        '-em', '--emmap_path',  help='Path to unsharpened EM map')
    locscale_emmap_input.add_argument(
        '-hm', '--halfmap_paths', help='Paths to first and second halfmaps', nargs=2)
    
    ## Input mask
    mask_input_group = parser.add_argument_group('Mask input arguments')
    mask_input_group.add_argument(
        '-ma', '--mask', help='Input filename mask')
    
    output_argument_group = parser.add_argument_group('Output arguments')
    output_argument_group.add_argument(
        '-v', '--verbose', help='Verbose output',action='store_true', default=False)
    output_argument_group.add_argument(
        '--print_report', help='Generate report PDF', action='store_true', default=False)
    output_argument_group.add_argument(
        '--report_filename', help='Filename for storing PDF output and statistics', default='locscale_report', type=str)
    output_argument_group.add_argument(
        '-op', '--output_processing_files', help='Path to store processing files', default=None, type=str)
    
    ## FDR parameters
    fdr_argument_group = parser.add_argument_group('FDR Confidence Mask arguments')
    fdr_argument_group.add_argument(
        '-fdr', '--fdr_threshold', help='FDR threshold for confidence mask', default=0.01, type=float)
    fdr_argument_group.add_argument(
        '-fdr_w', '--fdr_window_size', help='window size in pixels for FDR thresholding', default=None, type=int)
    fdr_argument_group.add_argument(
        '-avg_filter', '--averaging_filter_size', help='window size in pixels for FDR thresholding', default=3, type=int)
    fdr_argument_group.add_argument(
        '-fdr_f', '--fdr_filter', help='Pre-filter for FDR thresholding', default=None, type=float)
    fdr_argument_group.add_argument(
        '-th', '--mask_threshold', help='Threshold used to calculate the number of atoms and to decide the \
                                    envelope for initial placement of pseudo-atoms', default=0.99, type=float)

    ## Prediction parameters
    prediction_argument_group = parser.add_argument_group('Prediction arguments')
    prediction_argument_group.add_argument(
        '-model_path', '--model_path', help='Path to a custom trained model', default=None, type=str)
    prediction_argument_group.add_argument(
        '--use_low_context_model', help='Use a network which is trained on low context data (Model-Based LocScale targets)', \
                                    action='store_true', default=False,)
    prediction_argument_group.add_argument(
        '-bs', '--batch_size', help='Batch size for EMMERNET', default=8, type=int)
    prediction_argument_group.add_argument(
        '-gpus', '--gpu_ids', help="numbers of the selected GPUs, format: '1 2 3 ... 5'", required=False, nargs='+')
    prediction_argument_group.add_argument(
        '-cube_size','--cube_size', help='Size of the input cube for EMMERNET', default=32, type=int)

    
def add_locscale_arguments(locscale_parser):    
    ## Input model map file (mrc file) or atomic model (pdb file)
    locscale_parser.add_argument(
        '-o', '--outfile', help='Output filename', default="locscale_output.mrc")
    reference_map_group = locscale_parser.add_argument_group('Reference map arguments')
    model_input_group = reference_map_group.add_mutually_exclusive_group(required=False)
    model_input_group.add_argument(
        '-mm', '--model_map', help='Path to model map file')
    model_input_group.add_argument(
        '-mc', '--model_coordinates', help='Path to PDB file', default=None)
    reference_map_group.add_argument(
        '-mres', '--model_resolution', help='Resolution limit for Model Map generation', type=float)
    reference_map_group.add_argument(
        '-sym', '--symmetry', help='Impose symmetry condition for output', default='C1', type=str)
    
    ## LocScale main function parameters
    scaling_argument_group = locscale_parser.add_argument_group('Scaling arguments')
    scaling_argument_group.add_argument(
        '-wn', '--window_size', help='window size in pixels', default=None, type=int)
    scaling_argument_group.add_argument(
        '-mpi', '--mpi', help='MPI version', action='store_true', default=False)
    scaling_argument_group.add_argument(
        '-np', '--number_processes', help='Number of processes to use', type=int, default=1)

    ## Refinement parameters
    refinement_argument_group = locscale_parser.add_argument_group('Refinement arguments')
    refinement_argument_group.add_argument(
        '-ref_it', '--refmac_iterations', help='For atomic model refinement: number of refmac iterations', default=10, type=int)
    refinement_argument_group.add_argument(
        '-res', '--ref_resolution', help='Resolution target for Refmac refinement', type=float)
    refinement_argument_group.add_argument(
        '-p', '--apix', help='pixel size in Angstrom', type=float)
    refinement_argument_group.add_argument(
        '--add_blur', help='Globally sharpen the target map for REFMAC refinement', default=20, type=int)
    refinement_argument_group.add_argument(
        '--refmac5_path', help='Path to refmac5 executable', default=None, type=str)
    refinement_argument_group.add_argument(
        '--cref_pickle', help='Path for Cref filter for the target map of bfactor refinement', default=None, type=str)
    refinement_argument_group.add_argument(
        '-cif_info','--cif_info', help='Path to provide restrain information for refining the atomic model', default=None, type=str)    

    ## Integrated pseudo-atomic model method parameters
    hybrid_parser = locscale_parser.add_argument_group('Hybrid LocScale arguments')
    hybrid_parser.add_argument(
        '--complete_model', help='Add pseudo-atoms to areas of the map which are not modelled', action='store_true')
    hybrid_parser.add_argument(
        '-avg_w', '--averaging_window', help='Window size for filtering the fdr difference map for integrated pseudo-model',\
                                    default=3, type=int)

    ## Pseudo-atomic model method parameters
    pseudo_atomic_parser = locscale_parser.add_argument_group('Pseudo-atomic model arguments')
    pseudo_atomic_parser.add_argument(
        '--build_using_pseudomodel', help='Add pseudo-atoms to the map', action='store_true', default=False)
    pseudo_atomic_parser.add_argument(
        '-pm', '--pseudomodel_method', help='For pseudo-atomic model: method', default='gradient')
    pseudo_atomic_parser.add_argument(
        '-pm_it', '--total_iterations', help='For pseudo-atomic model: total iterations', default=50, type=int)
    pseudo_atomic_parser.add_argument(
        '-dst', '--distance', help='For pseudo-atomic model: typical distance between atoms', default=1.2, type=float)
    pseudo_atomic_parser.add_argument(
        '-mw', '--molecular_weight', help='Input molecular weight (in kDa)', default=None, type=float)
    pseudo_atomic_parser.add_argument(
        '--build_ca_only', help='For gradient pseudomodel building: use only Ca atoms with interatomic distance 3.8',\
                                    action='store_true',default=False)
    pseudo_atomic_parser.add_argument(
        '-s', '--smooth_factor', help='Smooth factor for merging profiles', default=0.3, type=float)
    pseudo_atomic_parser.add_argument(
        '--boost_secondary_structure', help='Amplify signal corresponding to secondary structures', default=1.5, type=float)
    
    ## Additional arguments for miscellaneous functions
    misc_parser = locscale_parser.add_argument_group('Miscellaneous arguments')
    misc_parser.add_argument(
        '--no_reference', help='Run locscale without using any reference information', action='store_true', default=False)
    misc_parser.add_argument(
        '--set_local_bfactor', help='For reference-less sharpening. Use this value to set the local b-factor of the maps',\
                                    type=float, default=20)
    misc_parser.add_argument(
        '--dev_mode', help='If true, this will force locscale to use the theoretical profile even if model map present \
                                    and will not check for user input consistency', action='store_true', default=False)
    misc_parser.add_argument(
        '--skip_refine', help='Ignore REFMAC refinement', action='store_true')


# **************************************************************************************
# ************************ Command line arguments EMMERNET *****************************
# **************************************************************************************

def add_emmernet_arguments(emmernet_parser):
    emmernet_parser.add_argument(
        '-o', '--outfile', help='Output filename', default="emmernet_output.mrc")
    emmernet_parser.add_argument(
        '-np', '--number_processes', help='Number of processes to use', type=int, default=1)
    emmernet_parser.add_argument(
        '-sym', '--symmetry', help='If not equal to C1, then symmetry averaging will be performed', default='C1', type=str)
    
    misc_parser = emmernet_parser.add_argument_group('Miscellaneous arguments')
    misc_parser.add_argument(
        '-no_mc','--no_monte_carlo', help='Disable Monte Carlo sampling of the output', action='store_true', default=False)
    misc_parser.add_argument(
        '-mc_it','--monte_carlo_iterations', help='Number of Monte Carlo iterations', default=15, type=int)
    misc_parser.add_argument(
        '-pb','--physics_based', help='Use physics-based model (under development!)', action='store_true')
    misc_parser.add_argument(
        '-download', '--download', help='Download the model weights', action='store_true', default=False)
    misc_parser.add_argument(
        '-s', '--stride', help='Stride for EMMERNET', default=16, type=int)
    

locscale_parser = argparse.ArgumentParser(prog="locscale",description="".join(description)) 
add_common_arguments(locscale_parser)
add_locscale_arguments(locscale_parser) 
## Add subparsers
sub_parser = locscale_parser.add_subparsers(dest='command')
## Add subparsers for feature_enhance
feature_enhance_parser = sub_parser.add_parser('feature_enhance', help='Enhance the features present in the input EM map through EMmerNet')
add_common_arguments(feature_enhance_parser)
add_emmernet_arguments(feature_enhance_parser)
version_parser = sub_parser.add_parser('version', help='Print version and exit')
test_parser = sub_parser.add_parser('test', help='Run all tests')
