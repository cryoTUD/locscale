# Script to automate LocScale for large number of files

from distutils.cmd import Command
from genericpath import isfile
import os
import sys

def get_defaults_dictionary():
    """
    Get the default input dictionary.
    """
    from locscale.main import locscale_parser
    defaults_dictionary = {}
    variables = locscale_parser._actions
    for variable in variables:
        defaults_dictionary[variable.dest] = variable.default
    
    return defaults_dictionary

## Create a class for LocScale inputs where the initial values are set to the defaults values using the get_defaults_dictionary() function
class LocScaleInputs:
    def __init__(self):
        self.input = get_defaults_dictionary()
        self.is_halfmap_input = None
        
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
        
        print("{} files copied".format(files_copied))

    def get_folder_from_input_paths(self):
        """
        Update the output folder from the input paths.
        """
        import os
        self.check_is_halfmap_input()

        if self.is_halfmap_input:
            halfmap_path_1 = os.path.abspath(self.input["halfmap_paths"][0])
            assert os.path.isfile(halfmap_path_1), "Halfmap 1 path is not valid"
            return os.path.dirname(halfmap_path_1)

        else:
            emmap_path = os.path.abspath(self.input["emmap_path"])
            assert os.path.isfile(emmap_path), "EM map path is not valid"
            return os.path.dirname(emmap_path)

        

class LocScaleRun:
    def __init__(self, Input, job_name, locscale_run_type, mpi_jobs, data_folder=None):
        self.input = Input
        self.job_name = job_name
        self.locscale_run_type = locscale_run_type
        self.command = ["locscale","run_locscale"]
        self.job_file_path = None
        self.job = None
        self.timeout = 4*3600 # 4 hours
        self.mpi_jobs = mpi_jobs
        self.input.check_is_halfmap_input()
        if data_folder is not None:
            self.data_folder = data_folder
        else:
            input_folder = self.input.get_folder_from_input_paths()
            self.data_folder = os.path.join(input_folder, job_name)
        
        if not os.path.exists(self.data_folder):
            os.mkdir(self.data_folder)
        
    def print_locscale_command(self):
        '''
        Print the command to the screen with a banner and run type. Command is a list of strings.
        '''
        print("="*80)
        print("Run type: {}".format(self.locscale_run_type))
        print("Command: \n")
        print(" ".join(self.command))
        print("="*80)

    def create_command(self):
        """
        Create locscale command from input dictionary
        """
        # Check for mandatory variables
        mandatory_variables_present = self.input.check_mandatory_variables(self.locscale_run_type)

        # Check for paths
        path_variables_valid = self.input.check_paths()
        print("Path variables valid: {}".format(path_variables_valid))
        print("Mandatory variables present: {}".format(mandatory_variables_present))
        
        if mandatory_variables_present and path_variables_valid:
            # Append unsharpened input based on is_halfmap_input
            if self.input.is_halfmap_input:
                self.command.append("--halfmap_paths")
                self.command.append(self.input.input["halfmap_paths"][0])
                self.command.append(self.input.input["halfmap_paths"][1])
            else:
                self.command.append("--emmap_path")
                self.command.append(self.input.input["emmap_path"])
            
            # For every other input, append to command list if value is not none or a boolean value 
            restricted_keys = ["halfmap_paths","emmap_path","model_coordinates","help",\
                "complete_model","dev_mode","averaging_window","mpi","symmetry","number_processes"]

            for key in self.input.input.keys():
                value = self.input.input[key]
                key_not_restricted = key not in restricted_keys
                value_is_not_none = value is not None
                value_is_boolean = isinstance(value, bool)
                if key_not_restricted and value_is_not_none and not value_is_boolean:
                    self.command.append("--{}".format(key))
                    self.command.append(str(value))

            # Append output folder

            self.command.append("--output_processing_files")
            processing_files_folder = os.path.join(self.data_folder, "processing_files_{}_{}".format(self.job_name, self.locscale_run_type))
            self.command.append(processing_files_folder)
            self.command.append("--verbose")

            if self.locscale_run_type == "model_based":
                self.command.append("--model_coordinates")
                self.command.append(self.input.input["model_coordinates"])
            
            if self.locscale_run_type == "model_free":
                self.command.append("--symmetry")
                self.command.append(str(self.input.input["symmetry"]))
            
            if self.locscale_run_type == "model_based_integrated":
                self.command.append("--complete_model")
                self.command.append("--model_coordinates")
                self.command.append(self.input.input["model_coordinates"])
                self.command.append("--averaging_window")
                self.command.append(str(self.input.input["averaging_window"]))
                self.command.append("--symmetry")
                self.command.append(str(self.input.input["symmetry"]))


                    
            
            # If MPI jobs are used, append the number of jobs to the command
            if self.mpi_jobs > 1:
                self.command.append("--mpi")
                self.command.insert(0,"mpirun")
                self.command.insert(1,"-np")
                self.command.insert(2,str(self.mpi_jobs))
            
                    
        else:
            print("Mandatory variables: {}".format(mandatory_variables_present))
            print("Path variables: {}".format(path_variables_valid))
            print("Command not created")
        
    
    def write_header_to_log_file(self,log_file_path):
        import os
        from datetime import datetime

        with open(log_file_path, "w") as f:
            f.write("="*80)
            f.write("\n")
            f.write("Run type: {}\n".format(self.locscale_run_type))
            f.write("User: {}\n".format(os.environ.get('USER')))
            f.write("Date: {}\n".format(datetime.now()))
            f.write("Command: \n")
            f.write(" ".join(self.command))
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
        self.create_command()
        self.write_header_to_log_file(log_file_path)
        self.print_locscale_command()

        # Create job file
        job = {
        "command": self.command, 
        "output_log": log_file_path, 
        "job_name": self.job_name,
        "timeout": int(self.timeout)}

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
        import os
        import subprocess

        self.fetch_job()

        log_file = open(self.job["output_log"], "a")
        # Submit job
        try:
            subprocess_output=subprocess.run(self.command, stdout=log_file, stderr=log_file, check=True, timeout=self.job["timeout"])
        except subprocess.CalledProcessError as exc:
            print("Error running LocScale for {}".format(self.job_name))
            print("Error: \n{}".format(exc))
            print("Command: \n{}".format(self.command))
            # Write error to log file
            log_file.write("Error running LocScale for {}\n".format(self.job_name))
            log_file.write("Error: \n{}\n".format(exc))
            log_file.close()
            return 2
        except subprocess.TimeoutExpired as exc:
            print("Timeout running LocScale for {}".format(self.job_name))
            print("Error: \n{}".format(exc))
            print("Time elapsed: {}".format(exc.timeout))
            print("Command: \n{}".format(self.command))
            # Write error to log file
            log_file.write("Timeout running LocScale for {}\n".format(self.job_name))
            log_file.write("Error: \n{}\n".format(exc))
            log_file.write("Time elapsed: {}\n".format(exc.timeout))
            log_file.close()
            return 1
        except Exception as exc:
            print("Error running LocScale for {}".format(self.job_name))
            print("Error: \n{}".format(exc))
            print("Command: \n{}".format(self.command))
            # Write error to log file
            log_file.write("Error running LocScale for {}\n".format(self.job_name))
            log_file.write("Error: \n{}\n".format(exc))
            log_file.close()
            return 420
        
        # Write success to log file
        log_file.write("Success running LocScale for {}\n".format(self.job_name))
        log_file.close()

        return 0
    
    def execute(self, dry_run=False):
        # print header
        self.print_header()
        # Prepare job
        self.prepare_job()
        # Submit job
        if not dry_run:
            returncode = self.submit_job()
        else:
            returncode = 10
            print("Dry run: job not submitted")
            self.print_footer()
            return returncode
        
        if returncode == 0:
            print("Successfully submitted job {}".format(self.job_name))
        elif returncode == 1:
            print("Timeout submitting job {}".format(self.job_name))
        elif returncode == 2:
            print("Error submitting job {}".format(self.job_name))
        else:
            print("Error submitting job {}".format(self.job_name))
        
        self.print_footer()
        return returncode

        
    
    def print_header(self):
        print("~"*80)
        print("Executing job name: {}".format(self.job_name))
        print("~"*80)
    
    def print_footer(self):
        print("~"*80)
        print("Finished job name: {}".format(self.job_name))
        print("~"*80)










            

        
            



            

        


        


        


        



        

