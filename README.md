[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/locscale)](https://pypi.org/project/locscale)
[![PyPI](https://img.shields.io/pypi/v/instamatic.svg?style=flat)](https://pypi.org/project/locscale/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/locscale)](https://pypi.org/project/locscale/)
[![DOI](https://zenodo.org/badge/DOI/10.7554/eLife.2713110.1007.svg)](https://doi.org/10.7554/eLife.27131)

# LocScale

`LocScale` is a program for reference-based density scaling (local sharpening) of cryo-EM maps with the aim to improve their interpretability.
  
`LocScale` is distributed as a portable stand-alone installation that includes all the needed libraries from: https://gitlab.tudelft.nl/aj-lab/locscale/releases.  

Please note that there is a GUI implemented version available as part of the [CCP-EM](http://www.ccpem.ac.uk/) project; it is also implemented in [Scipion](http://scipion.i2pc.es/).
<br>   

## Installation

If you use conda, create a new environment:

```
conda create -n locscale python=3.8 
conda activate locscale
```
0) Download the git repo: 
```
git clone https://gitlab.tudelft.nl/aj-lab/locscale.git
cd /path/to/repo/
```
1) Install packages using pip:

The setup.py file contains the list of packages and their versions which are used in this program. Use pip version 21.3 to ensure all packages and their version requirements are met. 

```bash
pip install -e . 
```
(this will create a pip module from your local git repository)

2) Install proshade and emda
```
pip install emda==1.1.3.post6
pip install proshade>=0.7.6.3
```
Note: 
(i) proshade and emda requires the following linux dependencies incase you have trouble installing the packages
(ii) You may need to update apt-get before installing the dependencies

```
sudo apt-get update
sudo apt-get upgrade
```

2a)
For proshade: (more information: https://github.com/michaltykac/proshade)

```
 sudo apt-get install gcc g++ make cmake git fftw3-dev liblapack-dev zlib1g-dev
```
2b) 
For emda:
```
sudo apt-get gfortran
```

3) Install GPU dependencies 
```
conda install -c anaconda cudatoolkit
conda install -c anaconda cudnn
```

4) Install MPI support
To run on multiple CPUs in parallel, `LocScale` uses MPI for distributed-memory parallelisation. Execution of the MPI-parallelized version requires installation of an MPI implementation (e.g. [openMPI](http://www.open-mpi.de/)) on your system.  

`LocScale` needs [mpi4py](http://pythonhosted.org/mpi4py/) to interface with the MPI implementation. Both can be installed via conda/pip:


```
conda install openmpi
pip install mpi4py
```

5) Run unittests 
(i) 
Change the active directory to "/path/to/locscale/tests/"
```
cd /path/to/locscale/tests/
```
(ii) LocScale tests
```
python -m unittest test_locscale.py -v
```
This should run the model-based and model-free locscale tests 

(iii) EMmerNet tests
```
python -m unittest test_emmernet.py -v
```
(iv) Symmetry tests (to check emda installation)
```
python -m unittest test_symmetry.py -v
```

Alternatively, download the portable installation with all libraries/dependencies included: https://gitlab.tudelft.nl/aj-lab/locscale/releases/latest.
<br> 

## How to use? 

Run LocScale (model-based) using the following syntax
```
python /path/to/locscale/locscale/main.py run_locscale -em path/to/emmap.mrc -mc path/to/model.pdb -res 3 -v -o model_based_locscale.mrc
```

Run LocScale (model-free) using the following syntax 
(this is the same syntax as above just not passing the model path runs the model free version automatically)
```
python /path/to/locscale/locscale/main.py run_locscale -em path/to/emmap.mrc -res 3 -v -o model_based_locscale.mrc
```

Run EMmerNet 
```
python /path/to/locscale/locscale/main.py run_emmernet -em path/to/emmap.mrc -v -trained_model model_based -gpus 0 -o emmernet_model_based.mrc
```

Note: For different EMmerNet models: Use the following syntax:
```
Model Based: -trained_model model_based
Model Free: -trained_model model_free
Ensemble Network: -trained_model ensemble
```

## Usage, tutorial and FAQs

Please see the [__Wiki__](https://gitlab.tudelft.nl/ajakobi/locscale/wikis/home) pages for usage instructions, FAQs and tutorial.
<br>  

## Credits

This project is using code from a number of third-party open-source projects. Projects used by `LocScale` are included under include/:

[FDRthresholding](https://git.embl.de/mbeckers/FDRthresholding) – tool for FDR-based density thresholding. License: 3-Clause BSD.  
[EMDA](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/emda/emda.html) – Python library for electron microscopy map and model manipulations. Licence: MPL-2    
[ProShade](https://github.com/michaltykac/proshade) - Protein Shape Description and Symmetry Detection. Licence: 3-Clause BSD

`LocScale` also makes use of [Refmac](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/content/refmac/refmac.html) – coordinate refinement program for macromolecular structures. Refmac is distributed as part of CCP-EM.

## Reference

If you found `LocScale` useful, please consider citing it:

- A. Bharadwaj and A.J. Jakobi, [Title still to be determined](do-not-know-yet), A Journal 1, 1-2 (2021)
- A.J. Jakobi, M. Wilmanns and C. Sachse, [Model-based local density sharpening of cryo-EM maps](https://doi.org/10.7554/eLife.27131), eLife 6: e27131 (2017).

---


For bug reports please use the [GitLab issue tracker](https://gitlab.tudelft.nl/aj-lab/locscale/issues).   
For questions or comments please contact <a.jakobi@tudelft.nl> or a.bharadwaj@tudelft.nl.
