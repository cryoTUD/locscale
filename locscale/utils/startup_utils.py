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

def launch_locscale_no_mpi(args):
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from locscale.utils.general import write_out_final_volume_window_back_if_required
    from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory
    import os 

    input_file_directory = get_input_file_directory(args) ## Get input file directory

    ## Print start
    start_time = datetime.now()
    print_start_banner(start_time, "LocScale")

    ## Check input
    check_user_input(args)   ## Check user inputs  
    if args.verbose:
        print_arguments(args)
    
    ## Change to output directory
    copied_args = change_directory(args, args.output_processing_files)  ## Copy the contents of files into a new directory
    parsed_inputs_dict = prepare_mask_and_maps_for_scaling(copied_args)
    ## Run LocScale non-MPI 
    LocScaleVol = run_window_function_including_scaling(parsed_inputs_dict)
    parsed_inputs_dict["output_directory"] = input_file_directory
    write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
    ## Print end
    print_end_banner(datetime.now(), start_time=start_time)


def launch_locscale_mpi(args):
    from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
    from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
    from locscale.utils.general import write_out_final_volume_window_back_if_required
    from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory
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
            parsed_inputs_dict["output_directory"] = input_file_directory
            write_out_final_volume_window_back_if_required(copied_args, LocScaleVol, parsed_inputs_dict)
            print_end_banner(datetime.now(), start_time=start_time)
    except Exception as e:
        print("Process {} failed with error: {}".format(rank, e))
        comm.Abort()
        raise e