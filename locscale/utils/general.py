import numpy as np
import math


def round_up_to_even(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if not return one higher
        return ceil_x
    else:
        return ceil_x+1

def round_up_to_odd(x):
    ceil_x = math.ceil(x)
    if ceil_x % 2 == 0:   ## check if it's even, if so return one higher
        return ceil_x+1
    else:
        return ceil_x

def true_percent_probability(n):
    x = np.random.uniform(low=0, high=100)
    if x <= n:
        return True
    else:
        return False

def copy_file_to_folder(full_path_to_file, new_folder):
    import shutil
    
    source = full_path_to_file
    file_name = source.split("/")[-1]
    destination = "/".join(new_folder.split("/")+[file_name])
    shutil.copyfile(source, destination)
    
    return destination

def change_directory(args, folder_name="processed"):
    import os    
    from locscale.utils.general import copy_file_to_folder
    
    current_directory = os.getcwd()
    new_directory = "/".join(current_directory.split("/")+[folder_name])
    os.mkdir(new_directory)
    
    for arg in vars(args):
        value = getattr(args, arg)
        if isinstance(value, str):
            if os.path.exists(value) and arg != "outfile":
                print("Copying {} to a new directory: \n{}".format(value, new_directory))
                new_location=copy_file_to_folder(value, new_directory)
                print("file saved at: {}".format(new_location))
                setattr(args, arg, new_location)
    
    os.chdir(new_directory)
    
    return args
                