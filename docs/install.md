# Installing LocScale 2.0 
!!! note "System Compatibility"
    LocScale should run on any CPU system with Linux, OS X or Windows subsytem for Linux (WSL). 
    <br><br>
    
    __GPU__: {==To run the `feature_enhance` option LocScale requires the availability of a GPU==}. It is possible to run it on CPUs but computation will be slow(er).  
    __OpenMPI__: LocScale allows parallelisation on multi-CPU environment with OpenMPI
    !!! warning "Installation on Apple silicon" 
        GPU support on Apple silicon (MX chip) is currently buggy. We are working on resolving this.  

We recommend to use [Conda](https://docs.conda.io/en/latest/) for a local working environment. See [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html#anaconda-or-miniconda) for more information on what Conda flavour may be the right choice for you, and [here](https://www.anaconda.com/products/distribution) for Conda installation instructions.

<div class="grid cards" markdown>

-   :material-clock-fast:{ .lg .middle } __Quick installation__

    ---

    Install [`LocScale 2.0`](#) and get up
    and running in minutes.
    
    [:octicons-arrow-right-24: Set up in 5 min](#quick)

-   :material-sticker-text-outline:{ .lg .middle } __Step-by-step instructions__

    ---

    Step-by-step installation installation instructions.

    [:octicons-arrow-right-24: Install guide](#detailed)
</div>


### Installation via environment files: {#quick}
##### 1. Install REFMAC5 via CCP4/CCPEM
LocScale needs a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure refmac5 is accessible to the program. 

##### 2. Install LocScale using environment files 

There are two yml files in the repo: environment_cpu.yml and environment_gpu.yml. We recommend you to download and install the GPU version.

Once you download the yml file of your choice: 
```bash
conda env create -f /path/to/environment_cpu.yml
conda activate cpu_locscale
```
or 
```bash
conda env create -f /path/to/environment_gpu.yml
conda activate gpu_locscale
```
### Installation using PyPi (pip) {#detailed}
You can also follow these steps to install locscale using pip.

##### 1. Create and activate a new conda environment

```bash title="1. Create and activate a new conda environment"
conda create -n locscale python=3.8 
conda activate locscale
```
##### 2. Install fortran compiler
LocScale uses Fortran code to perform symmetry operations and requires a Fortran compiler to be present in your system. You can install `gfortran` from conda-forge.
```bash
conda install -c conda-forge gfortran
```
##### 3. Install REFMAC5 via CCP4/CCPEM

The model-based and hybrid map sharpening modes of LocScale need a working instance of [REFMAC5](https://www2.mrc-lmb.cam.ac.uk/groups/murshudov/index.html). If you already have CCP4/CCPEM installed check if the path to run `refmac5` is present in your environment. For model-free sharpening and confidence-aware density modification REFMAC5 is not required. 

```bash
which refmac5
```

If no valid path is returned, please install [CCP4](https://www.ccp4.ac.uk/download/) to ensure REFMAC5 is accessible. 

##### 4. Install LocScale and dependencies using pip:

###### Recommended installation
We recommend using pip for installation. Use pip version 21.3 or later to ensure all packages and their version requirements are met. 

```bash
pip install locscale 
```

###### Install development version
If you would like to install the latest development version of locscale, use the following command to install from the git repository. 
```bash
pip install git+https://gitlab.tudelft.nl/aj-lab/locscale.git
```

To install the git repository in editable mode, clone the repository, navigate to the `locscale` directory, and run `pip install -e .`

##### 5. Testing

To test functionality after installation, you can run LocScale unit tests using the following command:

```bash
locscale test
```

### LocScale 2.0 in CCPEM Doppio
