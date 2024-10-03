#
# Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the program 'LocScale'
# written by the Author(s).
# Copyright (C) 2021 Alok Bharadwaj and Arjen J. Jakobi
# This software may be modified and distributed under the terms of the BSD license. 
# You should have received a copy of the BSD 3-clause license along with this program (see LICENSE file file for details).
# If not see https://opensource.org/license/bsd-3-clause/.
#

from __future__ import division, absolute_import, print_function
from setuptools.command.install import install
from setuptools.command.develop import develop
from numpy.distutils.core import Extension, setup
import pathlib
import os
import sys 

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))


def get_version():
    import pathlib    
    locscale_path = pathlib.Path(__file__).parent.resolve()
    version_text = (locscale_path / "locscale" / "__version__.py").read_text()
    version = version_text.split("=")[1][1:-1]
    return version
                      

def add_current_environment_to_init():
    import pathlib
    import os 
    locscale_path = pathlib.Path(__file__).parent.resolve()
    init_path = os.path.join(locscale_path, "locscale", "__init__.py")
    colab_env_str = str(os.environ.get('LOCSCALE_COLAB_ENV'))
    with open(init_path, "a") as f:
        f.write(f'\n__LOCSCALE_COLAB_ENV__ = {colab_env_str}\n')
        
def locscale_test_suite():
    import unittest
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='test_*.py')
    return test_suite

def check_for_refmac():
    import os
    from shutil import which
    
    refmac5_path = which("refmac5")
    if refmac5_path is None:
        print("Refmac5 is not installed. Some programs might not work.")
    else:
        print("Refmac5 is installed at {}".format(refmac5_path))
        print("If you want to use a different binary please use the --refmac5_path option or alias it to refmac5")


locscale_path = pathlib.Path(__file__).parent.resolve()
long_description = (locscale_path / "README.md").read_text()

class PostDevelopCommand(develop):
    """ Post-installation for development mode. """

    def run(self):
        develop.run(self)
        # Check if refmac5 is installed
        check_for_refmac()

class PostInstallCommand(install):
    """Post-installation for installation mode."""

    def run(self):
        install.run(self)
        # Check if refmac5 is installed
        check_for_refmac()
            
## Modify installation structure based on environment variables for different platforms
# LOCSCALE_COLAB_ENV
if os.getenv('LOCSCALE_COLAB_ENV'):
    ext_modules = []
else: 
    ex1 = Extension(name='fcodes_fast',
                    sources=['locscale/include/symmetry_emda/fcodes_fast.f90'])
    ext_modules = [ex1]

setup(name='locscale',
    version="2.2.5",
    author='Alok Bharadwaj, Arjen J. Jakobi, Reinier de Bruin',
    url='https://gitlab.tudelft.nl/aj-lab/locscale',
    project_urls={
        "Bug Tracker": "https://gitlab.tudelft.nl/aj-lab/locscale/issues",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: BSD License",
    ],
    description='Contrast optimization for cryo-EM maps',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='3-clause BSD',
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'locscale = locscale.main:main',
        ],
    },
    ext_modules=ext_modules,
    cmdclass={'develop': PostDevelopCommand,
                'install': PostInstallCommand},

    zip_safe=False)
