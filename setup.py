from __future__ import division, absolute_import, print_function
import setuptools
from setuptools.command.install import install
from setuptools.command.develop import develop

from numpy.distutils.core import setup, Extension

def compile_fcodes_fast():
  from numpy import f2py
  from locscale.utils.file_tools import get_locscale_path
  import os
  import subprocess

  fcodes_path = os.path.join(get_locscale_path(), "locscale","include", "fcodes_fast.f90")
  target_dir = os.path.dirname(fcodes_path)
  current_dir = os.getcwd()
  os.chdir(target_dir)
  subprocess.run(["f2py3.8","-c","-m","fcodes_fast",fcodes_path,"--build-dir",target_dir,"--quiet"])
  os.chdir(current_dir)
  

def download_emmernet_models():
  from locscale.utils.file_tools import get_locscale_path, download_emmernet_model_from_url, extract_tar_files_in_folder
  import os

  ## Create folder to download emmernet models
  emmernet_models_path = os.path.join(get_locscale_path(), "locscale","emmernet", "emmernet_models")
  if not os.path.exists(emmernet_models_path):
    os.makedirs(emmernet_models_path, exist_ok=True)
    download_emmernet_model_from_url(emmernet_models_path)
    extract_tar_files_in_folder(emmernet_models_path, use_same_folder=True)

def download_test_data():
  from locscale.utils.file_tools import get_locscale_path, download_test_data_from_url, extract_tar_files_in_folder
  import os

  ## Create folder to download tests_data
  test_data_path = os.path.join(get_locscale_path(), "tests","test_data")
  if not os.path.exists(test_data_path):
    os.makedirs(test_data_path, exist_ok=True)
    download_test_data_from_url(test_data_path)
    extract_tar_files_in_folder(test_data_path, use_same_folder=True)

def update_conda_environment():
  import subprocess
  from locscale.utils.file_tools import get_locscale_path
  import os
  subprocess.run(["conda", "env", "update", "--file", os.path.join(get_locscale_path(), "environment.yml")])

    
class PostDevelopCommand(develop):
  """ Post-installation for development mode. """
  def run(self):
    import os
    develop.run(self)   

    # Download emmernet models
    download_emmernet_models()

    # Download test data
    download_test_data()

    # Update conda environment
    update_conda_environment()

    # Compile fcodes_fast
    compile_fcodes_fast()

    
    

setup(name='locscale',
    version='2.0',
    author='Alok Bharadwaj, Arjen J. Jakobi, Reinier de Bruin and Carsten Sachse',
    url='https://gitlab.tudelft.nl/aj-lab/locscale',
    description= 'Contrast optimization for cryo-EM maps',
    license='3-clause BSD',
    packages=setuptools.find_packages(),
    install_requires=['matplotlib>=3.3.4','biopython>=1.78','numpy>=1.19.5','scipy>=1.5.4','pandas>=1.1.5',\
                      'mrcfile>=1.3.0','gemmi>=0.4.8','pypdb>=2.0','sklearn>=0.0','pwlf>=2.0.4','tqdm>=4.62.3',\
                      'more_itertools>=8.10.0','servalcat>=0.2.23','tensorflow==2.6.0','tensorflow-addons==0.14.0',\
                        'keras==2.6.0','tensorflow_datasets==4.5.2','pyfiglet>=0.8.post1'],
    entry_points={
      'console_scripts': [
          'locscale = locscale.main:main',
                          ],
      },
      cmdclass={'develop': PostDevelopCommand},
    zip_safe= False)

