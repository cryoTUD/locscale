import setuptools
from setuptools import setup
#from numpy.distutils.core import setup, Extension
import locscale

setup(name='locscale',
    version=locscale.__version__,
    author='Alok Bharadwaj and Arjen J. Jakobi',
    url='https://gitlab.tudelft.nl/aj-lab/locscale',
    description= 'Reference-based density scaling (sharpening) of cryo-EM maps',
    license='3-clause BSD',
    packages=setuptools.find_packages(),
    install_requires=['numpy>=1.19.5','scipy>=1.5.4','pandas>=1.1.5','mrcfile>=1.3.0','gemmi>=0.4.8','pypdb>=2.0','sklearn>=0.0','pwlf>=2.0.4','tqdm>=4.62.3','emda>=1.1.3','more_itertools>=8.10.0'],
    entry_points={
      'console_scripts': [
          'locscale = locscale.main:main',
                          ],
      },
    zip_safe= False)

