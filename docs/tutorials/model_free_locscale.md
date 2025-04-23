# Model-free LocScale<br><sup>Reference-based local sharpening without atomic model</sup>

If no atomic model is available, or if you do not want to use prior model information, you can use the model-free mode of LocScale. This method will predict a reference map using the ```EMmerNet``` network by default and is the recommended procedure for model-free local sharpening.  

<br>

!!!info "Model-free LocScale workflow"
    <br>
    ![alt text](img/model_free.png)

### Usage

```bash
locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -v -o model_free_locscale.mrc
```

Here, ```halfmap1.mrc``` and ```halfmap2.mrc``` should be the unsharpened and unfiltered half maps from yourr 3D refinement. If you wish to use the full map instead, use the following command:

```bash
locscale -em path/to/fullmap.mrc -mc path/to/model.pdb -v -o model_free_locscale.mrc
```

!!! note "Point group symmetry"
    If your map has point group symmetry, you need to specify the symmetry to force a symmetrised reference map for scaling. You can do
    this by specifying the required point group symmetry using the `-sym/--symmetry` flag, e.g. for D2:

    ```bash
    locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -v -sym D2 -o model_free_locscale.mrc
    ```

The output will be a locally sharpened map scaled according to the scale factors derived from the ```EMmerNet```-predicted reference map.

!!! warning "Recommended use of unfiltered input maps"
    Note that using unfiltered maps as input is essential. If using previously filtered maps, information beyond the spatial filter cutoff cannot be recovered.   

#### Model-free local sharpening with pseudomodels
Another option for model-free sharpening is use a full pseudotatom model. This can be enabled by passing the ```--build_using_pseudomodel``` flag when invoking ```LocScale```. This mode will estimate the molecular volume using statistical thresholding and generate a pseudo-atomic model within the thresholded boundary to approximate the distribution of atomic scatterers and estimate the local B-factor. It will then generate an average reference profile for local sharpening based on the experimental data and expected properties for electron scattering of biological macromolecules [2]. Use this if the default ```EMmerNet```-based reference map generation does not work well for your data (e.g. if the map is too noisy or if the map has very low resolution).

```bash
locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -v -o model_free_locscale.mrc --build_using_pseudomodel
```

Usually all default parameters for pseudomodel and reference profile generation are fine and we do not recommend to change them, but you can [change](https://gitlab.tudelft.nl/aj-lab/locscale/-/wikis/home/) them if you deem fit.

!!! note "Point group symmetry"
    If your map has point group symmetry, you need to specify the symmetry to force the pseudomodel generator for produce a symmetrised
    reference map for scaling. You can do this by specifying the required point group symmetry using the `-sym/--symmetry` flag, e.g.
    for D2:

    ```bash
    locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -v -sym D2 -o model_free_locscale.mrc
    ```

The output will be a locally sharpened map scaled according to the local scale factors derived from the ADP distribution of the hybrid pseudoatom model.

!!! tip "Speed-up computation on multiple CPUs"
    To speed up computation, you can use multiple CPUs if available. LocScale uses [OpenMPI](https://www.open-mpi.org/)/[`mpi4py`](https://mpi4py.readthedocs.io/en/stable/) for parallelisation, which should have been automatically set up during installation. You
    can run it as follows:

    ```bash
    mpirun -np 4 locscale -hm path/to/halfmap1.mrc path/to/halfmap2.mrc -v -o model_free_locscale.mrc -mpi
    ```
