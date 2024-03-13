# Script to automate LocScale for large number of files

import argparse
from distutils.cmd import Command
from genericpath import isfile
import os
from socket import timeout
import sys
from subprocess import PIPE
from locscale.utils.startup_utils import launch_locscale_no_mpi
import json
import pickle


def get_defaults_dictionary(program="locscale"):
    """
    Get the default input dictionary.
    """
    from locscale.utils.parse_utils import locscale_parser, feature_enhance_parser
    parser_to_choose = locscale_parser if program == "locscale" else feature_enhance_parser
    
    defaults_dictionary = {}
    variables = parser_to_choose._actions
    for variable in variables:
        # if variable.dest is not help then add to dictionary
        if variable.dest != 'help':
            defaults_dictionary[variable.dest] = variable.default
    
    return defaults_dictionary

## Create a class for LocScale inputs where the initial values are set to the defaults values using the get_defaults_dictionary() function
class LocScaleInputs:
    def __init__(self):
        self.input = get_defaults_dictionary()
        self.is_halfmap_input = None
        self.args = argparse.Namespace()
    
    def clone(self):
        """
        Clone the input dictionary.
        """
        clonedInput = LocScaleInputs()
        clonedInput.input = self.input.copy()
        clonedInput.is_halfmap_input = self.is_halfmap_input
        clonedInput.args = self.args
        return clonedInput
        
    def check_if_key_is_path(self, key, return_value=False):
        """
        Check if a key is a path.
        """
        if key != "halfmap_paths":
            value = self.input[key]
            value_is_not_none = value is not None
            value_is_str = isinstance(value, str)
            if value_is_not_none and value_is_str:
                if os.path.isfile(value):
                    return True
            return False
        else:
            return True

    def check_is_halfmap_input(self):
        """
        Check if the input is halfmap input.
        """
        halfmap_paths_present = self.input["halfmap_paths"] is not None
        self.is_halfmap_input = halfmap_paths_present
    def check_mandatory_variables(self, locscale_run_type):
        """
        Check if all mandatory variables are set.
        """
        # Check if input files are given

        # unsharpened maps inputs
        halfmap_paths_present = self.input["halfmap_paths"] is not None
        emmap_path_present = self.input["emmap_path"] is not None
        
        self.is_halfmap_input = halfmap_paths_present
        
        if halfmap_paths_present or emmap_path_present:
            unsharpened_input_present = True
        else:
            unsharpened_input_present = False

        resolution_present = self.input["ref_resolution"] is not None
        input_model_present = self.input["model_coordinates"] is not None
        complete_model = self.input["complete_model"] 

        if locscale_run_type == "model_based":
            mandatory_variables_present = resolution_present and input_model_present and unsharpened_input_present
        elif locscale_run_type == "model_free":
            mandatory_variables_present = resolution_present and unsharpened_input_present
        elif locscale_run_type == "model_based_integrated":
            mandatory_variables_present = resolution_present and input_model_present and unsharpened_input_present and complete_model
        else:
            print("LocScale run type not recognized")
        
        return mandatory_variables_present
    
    def check_paths(self):
        path_variables = ["halfmap_paths", "emmap_path", "model_coordinates", "mask"]
        path_variables_valid = True
        
        self.check_is_halfmap_input()
        if self.is_halfmap_input:
            halfmap_1_path = self.input["halfmap_paths"][0]
            halfmap_2_path = self.input["halfmap_paths"][1]
            halfmaps_valid = os.path.isfile(halfmap_1_path) and os.path.isfile(halfmap_2_path)
            path_variables_valid = halfmaps_valid and path_variables_valid
        else:
            emmap_path = self.input["emmap_path"]
            emmap_valid = os.path.isfile(emmap_path)
            path_variables_valid = emmap_valid and path_variables_valid
        if self.input["model_coordinates"] is not None:
            model_path = self.input["model_coordinates"]
            model_valid = os.path.isfile(model_path)
            path_variables_valid = model_valid and path_variables_valid
        if self.input["mask"] is not None:
            mask_path = self.input["mask"]
            mask_valid = os.path.isfile(mask_path)
            path_variables_valid = mask_valid and path_variables_valid

        return path_variables_valid
    
    def copy_files_to_new_folder(self, output_dir):
        """
        Copy files in a dictionary to the output directory.
        """
        import shutil
        files_copied = 0
        for key in self.input.keys():
            if key == "halfmap_paths":
                halfmap_path_1 = self.input["halfmap_paths"][0]
                halfmap_path_2 = self.input["halfmap_paths"][1]
                old_file_path_1 = os.path.abspath(halfmap_path_1)
                old_file_path_2 = os.path.abspath(halfmap_path_2)
                new_file_path_1 = os.path.join(output_dir, os.path.basename(old_file_path_1))
                new_file_path_2 = os.path.join(output_dir, os.path.basename(old_file_path_2))
                shutil.copy(old_file_path_1, new_file_path_1)
                shutil.copy(old_file_path_2, new_file_path_2)
                self.input["halfmap_paths"] = [new_file_path_1, new_file_path_2]
                files_copied += 2
            elif self.check_if_key_is_path(key):
                value = self.input[key]
                old_file_path = os.path.abspath(value)
                new_file_path = os.path.join(output_dir, os.path.basename(value))
                shutil.copy(old_file_path, new_file_path)
                self.input[key] = new_file_path
                files_copied += 1
            else:
                continue
        

    def get_folder_from_input_paths(self):
        """
        Update the output folder from the input paths.
        """
        import os
        self.check_is_halfmap_input()

        if self.is_halfmap_input:
            halfmap_path_1 = self.input["halfmap_paths"][0]
            assert os.path.isabs(halfmap_path_1)
            assert os.path.isfile(halfmap_path_1), "Halfmap 1 path is not absolute"
            return os.path.dirname(halfmap_path_1)

        else:
            emmap_path = os.path.abspath(self.input["emmap_path"])
            assert os.path.isfile(emmap_path), "EM map path is not valid"
            return os.path.dirname(emmap_path)

    def update_args(self):
        self.args.__dict__.update(self.input)
        

class LocScaleRun:
    def __init__(self, Input, job_name, locscale_run_type, mpi_jobs, data_folder=None):
        self.input = Input.clone()
    
        self.job_name = job_name
        self.locscale_run_type = locscale_run_type
        self.job_file_path = None
        self.job = None
        self.timeout = 8*3600
        self.mpi_jobs = mpi_jobs
        self.input.check_is_halfmap_input()
        if data_folder is not None:
            self.data_folder = data_folder
        else:
            input_folder = self.input.get_folder_from_input_paths()
            self.data_folder = os.path.join(input_folder, job_name)
        
        
        if not os.path.exists(self.data_folder):
            os.mkdir(self.data_folder)

    def write_header_to_log_file(self,log_file_path):
        import os
        from datetime import datetime

        with open(log_file_path, "w") as f:
            f.write("="*80)
            f.write("\n")
            f.write("Run type: {}\n".format(self.locscale_run_type))
            f.write("User: {}\n".format(os.environ.get('USER')))
            f.write("Date: {}\n".format(datetime.now()))
            f.write("Arguments: \n")
            f.write(self.input.args.__str__())
            f.write("\n")
            f.write("="*80)
        
    def prepare_job(self):
        import json

        # Create output folder for this job
        job_folder = self.data_folder
        
        # Copy files to new folder
        self.input.copy_files_to_new_folder(job_folder)

        # Make a locscale output log file
        log_file_path = os.path.join(job_folder, "locscale_output.log")

        # Create command
        self.input.update_args()
        self.write_header_to_log_file(log_file_path)

        # Create job file
        job = {
        "args": self.input.args.__dict__, 
        "output_log": log_file_path, 
        "job_name": self.job_name,
        "timeout": int(self.timeout),
        "data_folder": self.data_folder,}
        

        # Write job file
        self.job_file_path = os.path.join(job_folder, self.job_name + "_job.json")
        with open(self.job_file_path, "w") as f:
            json.dump(job, f)

        self.job = job
    
    def fetch_job(self):
        import json
        if self.job is None:
            with open(self.job_file_path, "r") as f:
                self.job = json.load(f)

    def submit_job(self):

        self.fetch_job()

        log_file = open(self.job["output_log"], "a")
        sys.stdout = log_file
        sys.stderr = log_file
        print("Submitting job {}".format(self.job_name))

        args = argparse.Namespace()
        args.__dict__.update(self.job["args"])

        # Run LocScale
        launch_locscale_no_mpi(args)

        print("Job {} submitted".format(self.job_name))
    

    def print_header(self):
        print("~"*80)
        print("Executing job name: {}".format(self.job_name))
        print("~"*80)
    
    def print_footer(self):
        print("~"*80)
        print("Finished job name: {}".format(self.job_name))
        print("~"*80)

class LocScaleMultiple:
    def __init__(self, input_jobs, num_process_scaling, num_process_preparation=1):
        self.input_jobs = input_jobs.copy()
        self.num_process_preparation = num_process_preparation
        self.num_process_scaling = num_process_scaling
        self.parsed_inputs_dictionary_list = {}
        self.locscale_volumes = {}
        self.dry_run = False
        self.split_jobs = False
        self.resume_scaling = False

        for job in self.input_jobs:
            job.prepare_job()
    
    def launch_locscale(self, job_id):
        from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory
        from locscale.utils.startup_utils import print_arguments, print_end_banner, print_start_banner, launch_locscale_no_mpi
        from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
        from datetime import datetime
        from contextlib import redirect_stdout

        ## Collect all preprocessing steps for all jobs
        args = argparse.Namespace()
        args.__dict__.update(self.input_jobs[job_id].job["args"])
        job_file_path = self.input_jobs[job_id].job_file_path
        log_file_path = self.input_jobs[job_id].job["output_log"]
        data_folder = self.input_jobs[job_id].job["data_folder"]
        parsed_dict_pickle = os.path.join(data_folder, "parsed_inputs.pickle")
        job_name = self.input_jobs[job_id].job["job_name"]
        with open(log_file_path, "a") as f:
            with redirect_stdout(f):
                if not self.dry_run:
                    try:
                        launch_locscale_no_mpi(args)
                    except Exception as e:
                        print("Error in job {}: {}".format(job_name, e))
                        with open(log_file_path, "a") as f:
                            f.write("Error in job {}: {}".format(job_name, e))
                            f.write("\n")
                        returncode = 1
                    
                    with open(log_file_path, "a") as f:
                        f.write("\n")
                        f.write("="*80)
                        f.write("\n")
                        f.write("Finished job {} at {}".format(job_name, datetime.now()))
                        f.write("\n")
                        f.write("="*80)
                        f.write("\n")
                else:
                    print("Dry run: skipping job {}".format(job_name))
                    with open(log_file_path, "a") as f:
                        f.write("\n")
                        f.write("="*80)
                        f.write("\n")
                        f.write("Finished job {} at {}".format(job_name, datetime.now()))
                        f.write("\n")
                        f.write("="*80)
                        f.write("\n")
        
        returncode = 0
        return returncode
       
        
    def get_parsed_inputs_dict(self, job_id):
        from locscale.utils.file_tools import change_directory, check_user_input, get_input_file_directory
        from locscale.utils.startup_utils import print_arguments, print_end_banner, print_start_banner
        from locscale.utils.prepare_inputs import prepare_mask_and_maps_for_scaling
        from datetime import datetime
        from contextlib import redirect_stdout

        ## Collect all preprocessing steps for all jobs
        args = argparse.Namespace()
        args.__dict__.update(self.input_jobs[job_id].job["args"])
        job_file_path = self.input_jobs[job_id].job_file_path
        log_file_path = self.input_jobs[job_id].job["output_log"]
        data_folder = self.input_jobs[job_id].job["data_folder"]
        parsed_dict_pickle = os.path.join(data_folder, "parsed_inputs.pickle")
        job_name = self.input_jobs[job_id].job["job_name"]
        
        # Print start banner
        print("Job ID: {} \t Job name: {} \t Begin preparing inputs".format(job_id, job_name))
        parsed_inputs_dict = {}
        parsed_inputs_dict["job_name"] = job_name
        parsed_inputs_dict["job_file_path"] = job_file_path
        parsed_inputs_dict["job_id"] = str(job_id)
        parsed_inputs_dict["input_file_directory"] = data_folder
        output_file_path = os.path.join(data_folder, args.outfile)
        parsed_inputs_dict["ouput_file_path"] = output_file_path

        with open(log_file_path, "w") as f:
            with redirect_stdout(f):
                
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
                
                if not self.dry_run:
                    try:
                        parsed_inputs_dict_from_args = prepare_mask_and_maps_for_scaling(copied_args)
                        parsed_inputs_dict.update(parsed_inputs_dict_from_args)
                        parsed_inputs_dict["parsing_success"] = str(True)
                    except Exception as e:
                        print("Error preparing inputs for job {}".format(job_id))
                        print(e)
                        parsed_inputs_dict["parsing_success"] = str(False)
                else:
                    print("Dry run: skipping preprocessing")              
                    parsed_inputs_dict["parsing_success"] = "dry_run"

                with open(parsed_dict_pickle, "wb") as f:
                        pickle.dump(parsed_inputs_dict, f)
                
                # dump parsed inputs dictionary to pickle file

                    

        print("Job ID: {} \t Job name: {} \t Finished preparing inputs".format(job_id, job_name))


    def scale_amplitudes(self, job_id):
        from locscale.utils.scaling_tools import run_window_function_including_scaling
        from locscale.include.emmer.ndimage.map_utils import save_as_mrc
        from contextlib import redirect_stdout
        
        job_name = self.input_jobs[job_id].job["job_name"]
        print("Job ID: {} \t Job name: {} \t Begin scaling".format(job_id, job_name))
        # print all keys from parsed_inputs_dictionary_list
        

        data_folder = self.input_jobs[job_id].job["data_folder"]
        parsed_dict_pickle = os.path.join(data_folder, "parsed_inputs.pickle")
        
        with open(parsed_dict_pickle, "rb") as f:
            parsed_inputs_dict = pickle.load(f)
        
       # print("SCALING: json file for job id {}: {}".format(job_id, parsed_dict_pickle))
       # print("Parsed inputs dict for job {}: {}".format(job_id, parsed_inputs_dict))
        
        job_file_path = self.input_jobs[job_id].job_file_path
        log_file_path = self.input_jobs[job_id].job["output_log"]

        with open(log_file_path, "a") as f:
            with redirect_stdout(f):
                print("Scaling job {}".format(job_id))
                print("="*80)
                if not self.dry_run:
                    if parsed_inputs_dict is not None:
                        try:
                            locscale_volume = run_window_function_including_scaling(parsed_inputs_dict)
                            save_as_mrc(locscale_volume, parsed_inputs_dict["ouput_file_path"], apix=parsed_inputs_dict["apix"])
                            returncode = 0
                            
                        except Exception as e:
                            print("Error scaling job {}".format(job_id))
                            print(e)
                            returncode = 2
                    else:
                        print("Problem preparing inputs for job {}".format(job_id))
                        returncode = 1
                else:
                    print("Dry run: skipping scaling")
                    returncode = 0
        print("Job ID: {} \t Job name: {} \t Finished scaling: {}".format(job_id, job_name, returncode))
        return returncode

    def launch(self):
        from locscale.utils.scaling_tools import run_window_function_including_scaling, run_window_function_including_scaling_mpi
        from locscale.utils.general import write_out_final_volume_window_back_if_required, merge_sequence_of_sequences, split_sequence_evenly
        import os 
        import joblib
        import numpy as np

        num_jobs = len(self.input_jobs)
        # Create a list of jobs to be submitted
        job_list = list(np.arange(num_jobs))

        # Submit jobs
        if self.split_jobs:

            # Run the preprocessing steps in parallel
            if self.resume_scaling:
                print("Resuming scaling from stored state")
            else:
                joblib.Parallel(n_jobs=self.num_process_preparation, backend="loky", timeout=6*3600)(
                    joblib.delayed(self.get_parsed_inputs_dict)(job_id) for job_id in job_list)

            ## Run the scaling in parallel
            result = joblib.Parallel(n_jobs=self.num_process_scaling, backend="loky",timeout=6*3600)(
                joblib.delayed(self.scale_amplitudes)(job_id) for job_id in job_list)
        else:
            result = joblib.Parallel(n_jobs=self.num_process_scaling, backend="loky",timeout=6*3600)(
                joblib.delayed(self.launch_contrast_enhance)(job_id) for job_id in job_list)

        
        # Print the results for each job
        for job_id in job_list:
            print("Job {} finished with status {}".format(job_id, result[job_id]))
            


            
        
        

    

    
        







            

        
            



            

        


        


        


        



        

