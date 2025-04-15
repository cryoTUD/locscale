# LocScale 2.0 <br><sup>Confidence-weighted cryoEM map optimisation</sup> 

`LocScale` is an automated program for local sharpening of cryo-EM maps with the aim to improve their interpretability. It utilises general properties inherent to electron scattering from biological macromolecules to restrain the sharpening filter. These can be provided either from an existing atomic model, or inferred directly from the experimental density map.

!!! note "What's new in LocScale 2.0"
    - Completely automated process for local map sharpening 
    - [Feature_enhance](#4-confidence-aware-density-modification): a confidence-aware density modification tool to enhance features in cryo-EM maps using the `EMmerNet` neural network.
    - [Hybrid sharpening](#2-run-locscale-using-a-partial-atomic-model): `LocScale` now supports reference-based sharpening when only partial atomic model information is available.
    - [Model-free sharpening](#3-run-locscale-without-atomic-model): `LocScale` now supports reference-based sharpening without the need to supply any atomic model information.
    - Full support for point group symmetry (helical symmetry to follow).

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

---

## What is LocScale2

The primary purpose of dlkkfklskflskfkle

What does the 'local scale'
