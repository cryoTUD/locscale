import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from datetime import datetime
import pyfiglet

progname = os.path.basename(sys.argv[0])
author = '\n\nAuthors: Arjen J. Jakobi (TU Delft), Alok Bharadwaj (TU Delft), Reinier de Bruin (TU Delft) \n\n'
version = progname + '  2.0'


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
locscale_parser.add_argument('--output_report', action='store_true', default=False,help='Print a PDF copy of the report')
locscale_parser.add_argument('--report_filename', type=str, help='Filename for storing PDF output and statistics', default="locscale_report")
locscale_parser.add_argument('-op', '--output_processing_files', type=str, help='Path to store processing files', default="processing_files")

## LocScale main function parameters
locscale_parser.add_argument('-wn', '--window_size', type=int, help='window size in pixels', default=None)
locscale_parser.add_argument('-mpi', '--mpi', action='store_true', default=False,help='MPI version')
locscale_parser.add_argument('-np', '--number_processes', help='Number of processes to use', type=int, default=1)

## Refinement parameters
locscale_parser.add_argument('-ref_it', '--refmac_iterations', help='For atomic model refinement: number of refmac iterations', default=15)
locscale_parser.add_argument('-res', '--ref_resolution', type=float, help='Resolution target for Refmac refinement')
locscale_parser.add_argument('-p', '--apix', type=float, help='pixel size in Angstrom')
locscale_parser.add_argument('--add_blur', type=int, help='Globally sharpen the target map for REFMAC refinement', default=20)
locscale_parser.add_argument('--refmac5_path', type=str, help='Path to refmac5 executable', default=None)

## Model map parameters
locscale_parser.add_argument('-mres', '--model_resolution', type=float, help='Resolution limit for Model Map generation')
locscale_parser.add_argument('-sym', '--symmetry', default='C1', type=str, help='Impose symmetry condition for output')

## FDR parameters
locscale_parser.add_argument('-fdr_w', '--fdr_window_size', type=int, help='window size in pixels for FDR thresholding', default=None)
locscale_parser.add_argument('-fdr_f', '--fdr_filter', type=float, help='Pre-filter for FDR thresholding', default=None)

## Pseudo-atomic model method parameters
locscale_parser.add_argument('-pm', '--pseudomodel_method', help='For pseudo-atomic model: method', default='gradient')
locscale_parser.add_argument('-pm_it', '--total_iterations', type=int, help='For pseudo-atomic model: total iterations', default=None)
locscale_parser.add_argument('-dst', '--distance', type=float, help='For pseudo-atomic model: typical distance between atoms', default=1.2)
locscale_parser.add_argument('-mw', '--molecular_weight', help='Input molecular weight (in kDa)', default=None)
locscale_parser.add_argument('--build_ca_only', help='For gradient pseudomodel building: use only Ca atoms with interatomic distance 3.8', action='store_true',default=False)
locscale_parser.add_argument('-s', '--smooth_factor', type=float, help='Smooth factor for merging profiles', default=0.3)
locscale_parser.add_argument('--boost_secondary_structure', type=float, help='Amplify signal corresponding to secondary structures', default=1.5)
locscale_parser.add_argument('--no_reference', action='store_true', default=False,help='Run locscale without using any reference information')

## non-default arguments
locscale_parser.add_argument('--dev_mode', action='store_true', default=False,help='If true, this will force locscale to use the theoretical profile even if model map present and will not check for user input consistency')
locscale_parser.add_argument('--ignore_profiles', help='Ignore average secondary structure profile during local scaling', action='store_true')
locscale_parser.add_argument('--skip_refine', help='Ignore REFMAC refinement', action='store_true')

# **************************************************************************************
# ************************ Command line arguments EMMERNET *******************************
# **************************************************************************************

## Input either unsharpened EM map or two halfmaps
emmernet_emmap_input = emmernet_parser.add_mutually_exclusive_group(required=True)
emmernet_emmap_input.add_argument('-em', '--emmap_path',  help='Path to unsharpened EM map')
emmernet_emmap_input.add_argument('-hm', '--halfmap_paths', help='Paths to first and second halfmaps')

## Output arguments
emmernet_parser.add_argument('-o', '--outfile', help='Output filename', default="emmernet_output.mrc")
emmernet_parser.add_argument('-op', '--output_processing_files', type=str, help='Path to store processing files', default="processing_files")
emmernet_parser.add_argument('-v', '--verbose', action='store_true',help='Verbose output')

## Emmernet main function parameters
emmernet_parser.add_argument('-trained_model','--trained_model', help='Type of emmernet model to use', \
                            choices=['model_based', 'model_free', 'ensemble'], default='model_based')
emmernet_parser.add_argument('-s', '--stride', help='Stride for EMMERNET', default=16, type=int)
emmernet_parser.add_argument('-bs', '--batch_size', type=int, help='Batch size for EMMERNET', default=8)
emmernet_parser.add_argument("-gpus", "--gpu_ids", nargs='+', help="numbers of the selected GPUs, format: '1 2 3 ... 5'", required=False)
emmernet_parser.add_argument('-download', '--download', help='Download the model weights', action='store_true', default=False)

############################################################################################
# ************************ Command line arguments TESTS ********************************** #
############################################################################################


def print_arguments(args):
    print("."*80)
    print('Input Arguments')
    print("."*80)
    for arg in vars(args):
        print('\t{}:  {}'.format(arg, getattr(args, arg)))        
    print("."*80)

def print_start_banner(start_time, text="Map Sharpening"):
    from textwrap import fill
    ## Definitions
    try:
        username = os.environ.get("USER")
    except:
        username = "Unknown"

    ## get today's date from start_time
    today_date = start_time.strftime("%d-%m-%Y")
    time_now = start_time.strftime("%H:%M:%S")

    ## Author credits
    
    if text == "LocScale":
        author_list = ["Arjen J. Jakobi (TU Delft)", "Alok Bharadwaj (TU Delft)"]
        contributor_list = ["Carsten Sachse (EMBL)"]
        version = "2.0"
    elif text == "EMmerNet":
        author_list = ["Arjen J. Jakobi (TU Delft)",  "Alok Bharadwaj (TU Delft)", "Reinier de Bruin (TU Delft)"]
        contributor_list = None
        version = "1.0"
    else:
        version = "x"

    ## Paper reference
    paper_ref_1 =  "Arjen J Jakobi, Matthias Wilmanns, Carsten Sachse (2017), \'Model-based local density sharpening of cryo-EM maps\', \'eLife 6:e27131\'"
    paper_ref_2 = "Alok Bharadwaj, Arjen J Jakobi (2022), \'Electron scattering properties of biological macromolecules and their use for cryo-EM map sharpening\', \'Faraday Discussions D2FD00078D\'"
    paper_ref_3 = "Alok Bharadwaj, Reinier de Bruin, Arjen J Jakobi (2022), \'TBD\'"
    print("="*80)
    print("="*80)
    result = pyfiglet.figlet_format(text, font = "big")
    print(result)
    print("\t"*6 + "Version: v{}".format(version))
    print("."*80)
    # Print user info and current time
    print("  |  ".join(["User: {}".format(username), "Date: {}".format(today_date), "Time: {}".format(time_now)]))
    print("\n")
    # Print author credits
    print("Authors:\n")
    for author in author_list:
        print("\t{} \n".format(author))
    # Print contributor credits if any
    if contributor_list is not None:
        print("Contributors:\n")
        for contributor in contributor_list:
            print("\t{} \n".format(contributor))
        
    # Print paper references
    print("References:\n")
    print(fill("{}".format(paper_ref_1), width=80, subsequent_indent="\t"))
    print(fill("{}".format(paper_ref_2), width=80, subsequent_indent="\t"))
    #print(wrap("{}".format(paper_ref_3), width=80))
    print("\n")
    if text == "EMmerNet":
        ## Print disclaimer for EMmerNet as this is in testing phase
        print("DISCLAIMER: Network Inpainting.\n")
        ## Print note on testing for network inpainting
        print(fill("EMmerNet is a neural network based map sharpening procedure. As such, there exists a risk of network hallucination " \
                +"i.e. the densities predicted by the network may not correspond to real densities. We are trying hard to mitigate "\
                +"this risk and we have undertaken a number of tests to ensure that network inpainting is not a problem. "\
                +"We have taken measures to ensure minimal bias exists in the training phase by using appropriate training targets."\
                +"If you encounter obvious problems, please report this to the authors. "+"\n"\
                +"Arjen Jakobi: a.jakobi@tudelft.nl  ", width=80))


    print("="*80)
    print("="*80)
    

def print_end_banner(time_now, start_time):
    print("."*80)
    ## print processing time in minutes
    print("Processing time: {:.2f} minutes".format((time_now-start_time).total_seconds()/60))
    print("="*80)
    print("Dank je wel!")
    print("="*80)
    

def launch_emmernet(args):
    from locscale.emmernet.utils import check_emmernet_inputs, check_and_save_output
    from locscale.utils.file_tools import change_directory
    from locscale.emmernet.prepare_inputs import prepare_inputs
    from locscale.emmernet.run_emmernet import run_emmernet
    
    ## Print start
    start_time = datetime.now()
    print_start_banner(start_time, "EMmerNet")
    if args.verbose:
        print_arguments(args)
    
    ## Check input
    check_emmernet_inputs(args)

    ## Change to output directory
    current_directory = os.getcwd()
    copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
    ## Prepare inputs
    input_dictionary = prepare_inputs(copied_args)
    ## Run EMMERNET
    emmernet_output_dictionary = run_emmernet(input_dictionary)
    ## Check and save output
    os.chdir(current_directory)
    check_and_save_output(input_dictionary, emmernet_output_dictionary)

    ## Print end
    print_end_banner(datetime.now(), start_time)


def launch_amplitude_scaling(args):
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from locscale.utils.general import write_out_final_volume_window_back_if_required
    from locscale.utils.file_tools import change_directory, check_user_input
    import os 

  
    current_directory = os.getcwd()
    
    if not args.mpi:
        ## Print start
        start_time = datetime.now()
        print_start_banner(start_time, "LocScale")

        ## Check input
        check_user_input(args)   ## Check user inputs    
        if args.verbose:
            print_arguments(args)
        
        ## Change to output directory
        copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
        ## Prepare inputs
        parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
        ## Run LocScale non-MPI 
        LocScaleVol = run_window_function_including_scaling(parsed_inputs_dict)
        ## Change to current directory and save output
        os.chdir(current_directory)
        write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
        ## Print end
        print_end_banner(datetime.now(), start_time=start_time)

    elif args.mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        ## If rank is 0, check and prepare inputs
        try:
            if rank==0:
                ## Print start
                start_time = datetime.now()
                print_start_banner(start_time, "LocScale")
                check_user_input(args)   ## Check user inputs
                if args.verbose:
                    print_arguments(args)
                    copied_args = change_directory(args, args.output_processing_files)
                parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
            else:
                parsed_inputs_dict = None
            
            ## Wait for inputs to be prepared by rank 0
            comm.barrier()
            ## Broadcast inputs to all ranks
            parsed_inputs_dict = comm.bcast(parsed_inputs_dict, root=0)           
            ## Run LocScale MPI
            LocScaleVol, rank = run_window_function_including_scaling_mpi(parsed_inputs_dict)
            ## Change to current directory and save output 
            if rank == 0:
                os.chdir(current_directory)
                write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
                print_end_banner(datetime.now(), start_time=start_time)
        except Exception as e:
            print(e)

def test_everything():
    from locscale.tests.utils import download_and_test_everything
    download_and_test_everything()

def main():
    main_args = main_parser.parse_args()
    launch_command = main_args.command

    if launch_command == 'run_emmernet':
        launch_emmernet(main_args)
    elif launch_command == 'run_locscale':
        launch_amplitude_scaling(main_args)
    elif launch_command == "test":
        test_everything()
    

if __name__ == '__main__':
    main()
