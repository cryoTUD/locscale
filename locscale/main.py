import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from datetime import datetime
#from mpi4py import MPI

progname = os.path.basename(sys.argv[0])
datmod = "2021-09-07"  # to be updated by gitlab after every commit
author = '\n\nAuthors: Arjen J. Jakobi, Carsten Sachse, EMBL and Alok Bharadwaj'
version = progname + '  0.2' + '  (;' + datmod+ ')'

simple_cmd = 'python main.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -w 25 -o loc_scaled.mrc'

description = "*** Computes contrast-enhanced cryo-EM maps by local amplitude scaling using a reference map *** \n \
Update: Now perform local amplitude scaling without an explicit reference to an atomic model"

cmdl_parser = argparse.ArgumentParser(
description='*** Computes contrast-enhanced cryo-EM maps by local amplitude scaling using a reference map.  ***\n \
Update: Now perform local amplitude scaling without an explicit reference to an atomic model' + \
('\nExample usage: \"{0}\". {1} on {2}'.format(simple_cmd, author, datmod)),formatter_class=RawTextHelpFormatter)


mpi_cmd = 'mpirun -np 4 python main.py -em emmap.mrc -mm modmap.mrc -ma mask.mrc -p 1.0 -w 10 -mpi -o scaled.mrc'

simple_cmd_use_pm = 'python main.py -em emmap.mrc -ma mask.mrc -o loc_scaled.mrc -use_pm'



cmdl_parser.add_argument('-em', '--em_map',  help='Input filename EM map')
cmdl_parser.add_argument('-hf1', '--half_map1', help='Input filename first half map')
cmdl_parser.add_argument('-hf2', '--half_map2', help='Input filename second half map')
cmdl_parser.add_argument('-mm', '--model_map', help='Input filename PDB map')
cmdl_parser.add_argument('-ma', '--mask', help='Input filename mask')
cmdl_parser.add_argument('-mc', '--model_coordinates', help='Input PDB files', default=None)
cmdl_parser.add_argument('-mw', '--molecular_weight', help='Input molecular weight (in kDa)', default=None)
cmdl_parser.add_argument('-o', '--outfile', help='Output filename')
cmdl_parser.add_argument('-op', '--output_processing_files', type=str, help='Path to store processing files', default="processing_files")
cmdl_parser.add_argument('-res', '--ref_resolution', type=float, help='Resolution target for Refmac refinement')
cmdl_parser.add_argument('-mres', '--model_resolution', type=float, help='Resolution limit for Model Map generation')
cmdl_parser.add_argument('-p', '--apix', type=float, help='pixel size in Angstrom')
cmdl_parser.add_argument('-wn', '--window_size', type=int, help='window size in pixels', default=None)
cmdl_parser.add_argument('-fdr_w', '--fdr_window_size', type=int, help='window size in pixels for FDR thresholding', default=None)
cmdl_parser.add_argument('-fdr_f', '--fdr_filter', type=float, help='Pre-filter for FDR thresholding', default=None)
cmdl_parser.add_argument('--ignore_profiles', help='Ignore average secondary structure profile during local scaling', action='store_true')
cmdl_parser.add_argument('--skip_refine', help='Ignore REFMAC refinement', action='store_true')
cmdl_parser.add_argument('-dst', '--distance', type=float, help='For pseudo-atomic model: typical distance between atoms', default=1.2)
cmdl_parser.add_argument('-pm', '--pseudomodel_method', help='For pseudo-atomic model: method', default='gradient')
cmdl_parser.add_argument('--build_ca_only', help='For gradient pseudomodel building: use only Ca atoms with interatomic distance 3.8', action='store_true',default=False)
cmdl_parser.add_argument('-pm_it', '--total_iterations', type=int, help='For pseudo-atomic model: total iterations', default=None)
cmdl_parser.add_argument('-ref_it', '--refmac_iterations', help='For pseudo-atomic model: number of refmac iterations', default=15)
cmdl_parser.add_argument('-sym', '--symmetry', default='C1', type=str, help='Impose symmetry condition for output')
cmdl_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,
                         help='MPI version call by: \"{0}\"'.format(mpi_cmd))
cmdl_parser.add_argument('--add_blur', type=int, help='Globally sharpen the map', default=0)
cmdl_parser.add_argument('-s', '--smooth_factor', type=float, help='Smooth factor for merging profiles', default=0.3)
cmdl_parser.add_argument('--boost_secondary_structure', type=float, help='Amplify signal corresponding to secondary structures', default=1.5)
cmdl_parser.add_argument('-v', '--verbose', action='store_true', default=False,
                         help='Verbose output')
cmdl_parser.add_argument('--dev_mode', action='store_true', default=False,
                         help='If true, this will force locscale to use the theoretical profile even if model map present and will not check for user input consistency')
cmdl_parser.add_argument('--report_filename', type=str, help='Filename for storing PDF output and statistics', default="locscale_report")
cmdl_parser.add_argument('--no_reference', action='store_true', default=False,
                         help='Run locscale without using any reference information')

def print_arguments(args):
    print('\n  LocScale Arguments\n')
    for arg in vars(args):
        print('    {} : {}'.format(arg, getattr(args, arg)))        
        
def print_start_banner(start_time):
    print("Launching amplitude scaling\n")
    print(start_time)
    print("**************************************************")
    print("*                      LOCSCALE                  *")
    print("**************************************************")

def print_end_banner(time_now, start_time):
    print("Finished amplitude scaling\n")
    print("Processing time: {}".format(time_now-start_time))
    print("**************************************************")
    print("*         Good luck with your research!          *")
    print("**************************************************")


def launch_amplitude_scaling(args):
    
    from locscale.utils.general import check_user_input
    
    
    #******************************************************************************************************************#
    
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi, write_out_final_volume_window_back_if_required
    from locscale.utils.general import change_directory
    import os 
    
    current_directory = os.getcwd()
    
    if not args.mpi:
        check_user_input(args)   ## Check user inputs
        start_time = datetime.now()
        print_start_banner(start_time)
        if args.verbose:
            print_arguments(args)
        copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
            
        parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
        
        LocScaleVol = run_window_function_including_scaling(parsed_inputs_dict)
        os.chdir(current_directory)
        write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
        
        '''
        ** INFORMATION ABOUT parsed_inputs_dict VARIABLE **
        parsed_inputs_dict.keys()
        ['emmap'] : numpy.ndarray   | Experimental Map
        ['modmap'] : numpy.ndarray  | Reference Map for scaling
        ['mask'] : numpy.ndarray    | Mask 
        ['wn'] :  int               | Window size, LocScale (pix) 
        ['apix'] : float            | Pixel size 
        ['use_theoretical'] : bool  | Flag to use theoretical       profiles
        ['scale_factor_args'] : dict| Parameters for scale factor computation
        
        ['verbose'] : bool          | Verbose parameter
        ['win_bleed_pad'] : bool    | Window bleed and pad         
        
        ['bfactor_info'] : list     | Information about bfactors, slopes and breakpoints from PWLF
        ['fsc_resolution'] : float  | FSC resolution
        ['PWLF_fit'] : float        | R^2 of PWLF for the radial profile of the input EM map
        ['emmap_path'] : str        | Em map path (after running mapmask command)
        ['mask_path'] : str         | Mask path (after running mapmask command)
        
        '''
        
        print("You can find the scaled map here: {}".format(args.outfile))
        print_end_banner(datetime.now(), start_time=start_time)
        
        

            
    
    elif args.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        
        try:
            if rank==0:
                check_user_input(args)   ## Check user inputs
                start_time = datetime.now()
                print_start_banner(start_time)
                if args.verbose:
                    print_arguments(args)
                    copied_args = change_directory(args, args.output_processing_files)
                
                parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
            
            else:
                parsed_inputs_dict = None
            
            comm.barrier()
            
            parsed_inputs_dict = comm.bcast(parsed_inputs_dict, root=0)
        
            #input_to_scaling = parsed_arguments[:-1]
            
            LocScaleVol, rank = run_window_function_including_scaling_mpi(parsed_inputs_dict)
            
            if rank == 0:
                os.chdir(current_directory)
                write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
        
                print("You can find the scaled map here: {}\n".format(args.outfile))
                print_end_banner(datetime.now(), start_time=start_time)
        except Exception as e:
            print(e)
            
        

def main():
    args = cmdl_parser.parse_args()

    launch_amplitude_scaling(args)

if __name__ == '__main__':
    main()
