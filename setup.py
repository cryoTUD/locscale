from __future__ import division, absolute_import, print_function
import setuptools
from numpy.distutils.core import setup, Extension


setup(name='locscale',
    version='2.0',
    author='Alok Bharadwaj, Arjen J. Jakobi, Reinier de Bruin and Carsten Sachse',
    url='https://gitlab.tudelft.nl/aj-lab/locscale',
    description= 'Contrast optimization for cryo-EM maps',
    license='3-clause BSD',
    packages=setuptools.find_packages(),
    install_requires=['matplotlib>=3.3.4','biopython>=1.78','numpy>=1.19.5','scipy>=1.5.4','pandas>=1.1.5','mrcfile>=1.3.0','gemmi>=0.4.8','pypdb>=2.0','sklearn>=0.0','pwlf>=2.0.4','tqdm>=4.62.3','more_itertools>=8.10.0','servalcat>=0.2.23','tensorflow==2.6.0','tensorflow-addons==0.14.0','keras==2.6.0','pyfiglet>=0.8.post1'],
    entry_points={
      'console_scripts': [
          'locscale = locscale.main:main',
                          ],
      },
    zip_safe= False)

