# LocScale 2.0 <br><sup>Confidence-weighted cryoEM map optimisation</sup> 
<div style="text-align: justify"><code>LocScale 2.0</code> is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.
</div>

!!! info "What's new in LocScale 2.0"  
    - Completely automated process for local map sharpening   
    - [Feature_enhance](#4-confidence-aware-density-modification): Confidence-weighted map optimisation by variational inference.
    - [Hybrid sharpening](#2-run-locscale-using-a-partial-atomic-model): `LocScale` now supports reference-based sharpening with partial (incomplete) models.  
    - [Model-free sharpening](#3-run-locscale-without-atomic-model): `LocScale` now supports reference-based sharpening without any model information.  
    - Full support for point group symmetry (helical symmetry to follow).  

---
## How does LocScale 2.0 work?

<div style="text-align: justify">
<code>LocScale 2.0</code> integrates physics-informed and deep learning-based map optimisation. Physical priors for map sharpening are based on established knowledge about expected properties of electron scattering by biological macromolecules <a href="#bottom">[1-3]</a>. Alternatively, a deep convolutional neural network [EMmerNet] trained on pairs of unmodified maps and maps optimised with the physics-informed scaling procedure can be used to predict local scale estimates. This information is then fed into a windowed scaling procedure to produce locally sharpened maps. Both of these workflows are map sharpening procedures operate in Fourier space, where structure factor amplitudes are corrected but map phases are locally unchanged. In a third workflow [<code>locscale_feature_enhance</code>], a Bayesian-approximate implementation of EMmerNet is used to predict an optimised map (which we call <b>feature-enhanced map</b>) along with its uncertainties. This procedure operates in real space and affects phases akin to density modification. LocScale 2.0 computes a voxel-wise confidence score that quantifies the uncertainty of this prediction, which can be mapped onto the map to guide interpretation. 
</div>
<br>
<div style="display: flex; flex-direction: column; align-items: center;">
  <div class="c-compare" style="--value:50%; position: relative; width: 600px; height: 400px; overflow: hidden;">
    <img class="c-compare__left"
         src="imgs/emd19995.png"
         alt="Raw map"
         style="position: absolute; top: 0; left: 0; width: 100%; height: 100%; object-fit: contain;" />

    <img class="c-compare__right"
         src="imgs/emd19995_fem.png"
         alt="Feature-enhanced map"
         style="position: absolute; top: 0; left: 0; width: 100%; height: 100%; object-fit: contain; clip-path: inset(0 0 0 var(--value));" />

    <input type="range" class="c-compare__range" min="0" max="100" value="50"
           oninput="this.parentNode.style.setProperty('--value', this.value + '%')"
           style="position: absolute; bottom: 10px; left: 10px; width: 90%; z-index: 10;" />
  </div>
</div>
<br>

!!! info "What are we looking at here?"
    This figure shows the effect of map optimisation with ```LocScale 2.0```, in this case using the ```feature_enhance``` option. ```LocScale 2.0``` attempts to simultaneously enhance high-resolution structural detail and improve contrast of low(er) resolution map regions associated with flexible domains, subunits with partial occupancy and contextual structures such as detergent micelles.

!!! warning 
    Confidence scores are not a 100% safeguard against misleading interpretations.

## Which map optimisation procedure should I use?


![alt](imgs/overview_methods.png)

In general we recommend using ```locscale_feature_enhance``` for map optimisation in ```LocScale 2.0``` whenever applicable. We have found this procedure to work robustly in a majority of cases and to provide the best compromise in preserving high-resolution detail and enhancing contrast of flexible or lower occupancy regions and contextual structures such as micelles.  


