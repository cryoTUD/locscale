import pathlib

locscale_path=pathlib.Path(__file__).parent.resolve()


def download_emmernet_model_from_url(download_folder):
    import wget
   
    url_model_based_emmernet = "https://surfdrive.surf.nl/files/index.php/s/HxRLgoZFYQEbf8Z/download"
    wget.download(url_model_based_emmernet, download_folder)

def download_test_data_from_url(download_folder):
    import wget
   
    #url_test_data = "https://surfdrive.surf.nl/files/index.php/s/xJKxGXR0LWGBDWM/download"
    url_test_data = "https://surfdrive.surf.nl/files/index.php/s/lk9CdNO5gszFll1/download"
    
    wget.download(url_test_data, download_folder)

def extract_tar_files_in_folder(tar_folder, use_same_folder=True):
    import tarfile
    import os
    if use_same_folder == 0:
        target_folder = tar_folder
    else:
        target_folder = os.path.dirname(tar_folder)

    for file in os.listdir(tar_folder):
        if file.endswith(".tar.gz"):
            print("Extracting: {}".format(file))
            tar = tarfile.open(os.path.join(tar_folder,file))
            tar.extractall(target_folder)
            tar.close()

def compile_fcodes_fast():
  from numpy import f2py
  import os
  import subprocess
  from shutil import which

  
  gfortran_location = which("gfortran")
  assert gfortran_location is not None, "gfortran is not installed"

  target_dir = os.path.dirname(fcodes_path)
  current_dir = os.getcwd()
  os.chdir(target_dir)
  subprocess.run(["f2py3","-c","-m","fcodes_fast",fcodes_path,"--quiet","--f90exec={}".format(gfortran_location)])
  os.chdir(current_dir)
  

def download_emmernet_models():
  import os
  
  ## Find the absolute path of the locscale folder
  import locscale
  locscale_path = os.path.dirname(locscale.__file__)
  print("locscale_path: {}".format(locscale_path))
  ## Create folder to download emmernet models
  emmernet_models_path = os.path.join(locscale_path, "locscale", "emmernet", "emmernet_models")
  if not os.path.exists(emmernet_models_path):
    os.makedirs(emmernet_models_path, exist_ok=True)
    download_emmernet_model_from_url(emmernet_models_path)
    extract_tar_files_in_folder(emmernet_models_path, use_same_folder=True)

def download_test_data():
  import os

  ## Create folder to download tests_data
  test_data_path = os.path.join(locscale_path, "tests","test_data")
  if not os.path.exists(test_data_path):
    os.makedirs(test_data_path, exist_ok=True)
    download_test_data_from_url(test_data_path)
    extract_tar_files_in_folder(test_data_path, use_same_folder=True)

def update_conda_environment():
  import subprocess
  import os
  # Install cudatoolkit using conda 
  subprocess.run(["conda", "install", "-c", "anaconda", "cudatoolkit","--yes"])

  # Install cudnn
  subprocess.run(["conda", "install", "-c", "anaconda", "cudnn","--yes"])

  # Install openmpi
  subprocess.run(["conda", "install", "-c", "conda-forge", "openmpi","--yes"])

  # Install mpi4py
  subprocess.run(["conda", "install", "-c", "conda-forge", "mpi4py==3.0.0","--yes"])

def check_refmac5_installed():
  from shutil import which
  # Check if refmac5 is installed
  refmac5_location = which("refmac5")
  if refmac5_location is None:
    raise UserWarning("Refmac5 is not installed. Please install it and try again.")
  else:
    print("Refmac5 is installed at {}".format(refmac5_location))
    print("If you want to use a different binary please use the --refmac5_location option")

def check_and_install_gfortran():
  from shutil import which
  from subprocess import run
  # Check if gfortran is installed
  gfortran_location = which("gfortran")
  if gfortran_location is None:
    print("Installing gfortran using conda")
    run(["conda", "install","-c","conda-forge","gfortran"])


def locscale_test_suite():
  test_loader = unittest.TestLoader()
  test_suite = test_loader.discover('tests', pattern='test_*.py')
  return test_suite

def run_locscale_tests():
  test_runner = unittest.TextTestRunner()
  test_suite = locscale_test_suite()
  test_runner.run(test_suite)

def main():

    # Download emmernet models
    download_emmernet_models()

    # Download test data
    download_test_data()

    # Update conda environment
    update_conda_environment()

    # Install gfortran if not installed
    check_and_install_gfortran()

    # Compile fcodes_fast
    # compile_fcodes_fast()

    # Check if refmac5 is installed
    check_refmac5_installed()

    # Run tests
    run_locscale_tests()

if __name__ == "__main__":
    main()
