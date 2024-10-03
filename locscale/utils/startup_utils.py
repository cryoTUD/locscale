#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license.  You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details). If not see https://opensource.org/license/bsd-3-clause/.
#

import os
import pyfiglet
from datetime import datetime
import sys
import locscale 

def print_arguments(args):
    print("."*80)
    print('Input Arguments')
    print("."*80)
    for arg in vars(args):
        print('\t{}:  {}'.format(arg, getattr(args, arg)))        
    print("."*80)

def print_start_banner(start_time, text="Map Sharpening"):
    from textwrap import fill
    import time
    ## Definitions
    try:
        username = os.environ.get("USER")
    except:
        username = "Unknown"

    disclaimer_text_path = os.path.join(os.path.dirname(locscale.__file__), "disclaimer.txt")
    disclaimer_text = open(disclaimer_text_path, "r").read()
                                        
    ## get today's date from start_time
    today_date = start_time.strftime("%d-%m-%Y")
    time_now = start_time.strftime("%H:%M:%S")

    ## Author credits
    
    if text == "LocScale":
        author_list = ["Arjen J. Jakobi (TU Delft)", "Alok Bharadwaj (TU Delft)"]
        contributor_list = ["Carsten Sachse (EMBL)"]
        version = locscale.__version__
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
        print(fill(disclaimer_text, width=80))
        # Sleep for 5 seconds
        time.sleep(5)


    print("="*80)
    print("="*80)
    

def print_end_banner(time_now, start_time):
    print("."*80)
    ## print processing time in minutes
    print("Processing time: {:.2f} minutes".format((time_now-start_time).total_seconds()/60))
    print("="*80)
    print("Dank je wel!")
    print("="*80)

def launch_locscale_no_mpi(args):
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from locscale.utils.general import write_out_final_volume_window_back_if_required
    from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory, set_modality_based_on_input, pretty_print_dictionary
    import os 

    input_file_directory = get_input_file_directory(args) ## Get input file directory

    ## Print start
    start_time = datetime.now()
    if args.verbose:
        print_start_banner(start_time, "LocScale")

    ## Check input
    check_user_input(args)   ## Check user inputs  
    args=set_modality_based_on_input(args) ## Set modality based on input
    if args.verbose:
        print_arguments(args)
    
    ## Change to output directory
    copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
    copied_args.logger.info("Running LocScale")
    copied_args.logger.info("-"*80)
    copied_args.logger.info(f"Arguments used: \n{pretty_print_dictionary(vars(copied_args))}")
    copied_args.logger.info("-"*80)
    parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
    ## Run LocScale non-MPI 
    LocScaleVol = run_window_function_including_scaling(parsed_inputs_dict)
    parsed_inputs_dict["output_directory"] = input_file_directory
    write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
    ## Print end
    if args.verbose:
        print_end_banner(datetime.now(), start_time=start_time)


def launch_locscale_mpi(args):
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from locscale.utils.general import write_out_final_volume_window_back_if_required
    from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory, set_modality_based_on_input, setup_logger
    import os 

    input_file_directory = get_input_file_directory(args) ## Get input file directory

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    ## If rank is 0, check and prepare inputs
    try:
        if rank==0:
            ## Print start
            start_time = datetime.now()
            if args.verbose:
                print_start_banner(start_time, "LocScale")
            check_user_input(args)   ## Check user inputs
            args=set_modality_based_on_input(args) ## Set modality based on input
            if args.verbose:
                print_arguments(args)
            copied_args = change_directory(args, args.output_processing_files)
            copied_args.logger.info("Running LocScale MPI")
            copied_args.logger.info("-"*80)
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
            parsed_inputs_dict["output_directory"] = input_file_directory
            write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
            if args.verbose:
                print_end_banner(datetime.now(), start_time=start_time)
    except Exception as e:
        print("Process {} failed with error: {}".format(rank, e))
        comm.Abort()
        raise e

def launch_contrast_enhance(args):
    from locscale.utils.startup_utils import launch_locscale_no_mpi, launch_locscale_mpi
    if args.mpi:
        launch_locscale_mpi(args)
    else:
        launch_locscale_no_mpi(args)

def launch_feature_enhance(args):
    from locscale.emmernet.utils import check_emmernet_inputs, check_and_save_output
    from locscale.utils.file_tools import change_directory, pretty_print_dictionary
    from locscale.utils.general import try_to
    from locscale.utils.plot_tools import compute_probability_distribution
    from locscale.emmernet.prepare_inputs import prepare_inputs
    from locscale.emmernet.run_emmernet import run_emmernet
    from locscale.emmernet.emmernet_functions import calculate_significance_map_from_emmernet_output
    
    ## Print start
    start_time = datetime.now()
    print_start_banner(start_time, "EMmerNet")
    if args.verbose:
        print_arguments(args)
    
    ## Check input
    check_emmernet_inputs(args)
    
    copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
    copied_args.logger.info("Running Feature Enhancement using EMmerNet")
    copied_args.logger.info("-"*80)
    copied_args.logger.info(f"Arguments used: \n{pretty_print_dictionary(vars(copied_args))}")
    copied_args.logger.info("-"*80)
    
    ## Prepare inputs
    input_dictionary = prepare_inputs(copied_args)
    
    ## Run EMMERNET
    emmernet_output_dictionary = run_emmernet(input_dictionary)
    emmernet_output_dictionary = check_and_save_output(input_dictionary, emmernet_output_dictionary)
    
    ## Run LocScale using EMmerNet output
    locscale_args = get_locscale_inputs_from_emmernet(input_dictionary, emmernet_output_dictionary)
    launch_contrast_enhance(locscale_args)
    
    # Calculate the p-values from the output of EMmerNet and LocScale
    locscale_output_filename = locscale_args.outfile 
    emmernet_output_mean_filename = emmernet_output_dictionary["output_filename_mean"]
    emmernet_output_var_filename = emmernet_output_dictionary["output_filename_var"]
    calculate_significance_map_from_emmernet_output(
        locscale_output_filename, emmernet_output_mean_filename, emmernet_output_var_filename, n_samples=input_dictionary["monte_carlo_iterations"])
    # Compute the probability distribution for calibration
    
    probability_args = {"locscale_path" : locscale_output_filename, 
                        "mean_prediction_path" : emmernet_output_mean_filename,
                        "var_prediction_path" : emmernet_output_var_filename,
                        "n_samples" :input_dictionary["monte_carlo_iterations"], 
                        "processing_files_folder" : input_dictionary["output_processing_files"]}
    
    try_to(compute_probability_distribution,**probability_args)
    
    print("EMmerNet finished successfully")
    print("Please check the pVDDT_scores.mrc to observe significantly hallucinated features")
    ## Print end
    print_end_banner(datetime.now(), start_time)
    
def get_locscale_inputs_from_emmernet(parsed_inputs, emmernet_output):
    from locscale.automate.tools import get_defaults_dictionary
    import argparse
    
    ## Get defaults dictionary for LocScale
    input_folder = parsed_inputs["input_folder"]
    output_filename = os.path.basename(emmernet_output["output_filename_for_locscale"])
    output_path = os.path.join(input_folder, output_filename)
    defaults_dictionary = get_defaults_dictionary("locscale")
    defaults_dictionary["emmap_path"] = parsed_inputs["xyz_emmap_path"]
    defaults_dictionary["mask"] = parsed_inputs["xyz_mask_path"]
    defaults_dictionary["model_map"] = emmernet_output["reference_map_for_locscale"]
    defaults_dictionary["verbose"] = False 
    defaults_dictionary["outfile"] = output_path
    defaults_dictionary["output_processing_files"] = parsed_inputs["output_processing_files"]
    defaults_dictionary["logger"] = parsed_inputs["logger"]
    defaults_dictionary["number_processes"] = parsed_inputs["number_processes"]
    
    locscale_args = argparse.Namespace(**defaults_dictionary)
    
    return locscale_args

def run_housekeeping():
    import sys 

    # Add installation date to __init__.py 
    add_installation_date()
    # Check if help message needs to be printed
    check_for_help_message(sys.argv)

def check_for_help_message(system_arguments):
    import sys 
    from locscale.utils.parse_utils import locscale_parser
    # Checks whether help is needed 
    if len(system_arguments) == 1:
        locscale_parser.print_help()
        sys.exit(1)
    
    else:
        launch_command = locscale_parser.parse_args().command
        if launch_command == 'feature_enhance':
            if len(system_arguments) == 2:
                locscale_parser.print_help()
                sys.exit(1)
        elif launch_command == 'version':
            pass
        elif launch_command == 'test':
            pass
        elif launch_command is None:
            pass
        else:
            raise ValueError("Unknown command: ", launch_command)
        
    
def add_installation_date():

    from datetime import datetime
    from locscale.utils.file_tools import get_locscale_path

    init_path = os.path.join(get_locscale_path(), "locscale","__init__.py")    
    # readlines
    with open(init_path, "r") as f:
        lines = f.readlines()

    # check if __installation_date__ is already present
    installation_date_added = False
    for line in lines:
        if "__installation_date__" in line:
            installation_date_added = True
            break
        
    # write __installation_date__ to __init__.py if not present
    if not installation_date_added:
        with open(init_path, "a") as f:
            f.write(f'\n__installation_date__ = "{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}"\n')
    else:
        pass 

        
def test_everything():
    from locscale.tests.utils import download_and_test_everything
    download_and_test_everything()

def print_version():
    run_housekeeping()
    print("LocScale")
    print("Version: ", locscale.__version__)
    try:
        print("Installed on: ", locscale.__installation_date__)
    except AttributeError:
        print("Installation date not available")

    print("Authors: Arjen J. Jakobi (TU Delft), Alok Bharadwaj (TU Delft), Reinier de Bruin (TU Delft)")
    print("Python version: ", sys.version)

