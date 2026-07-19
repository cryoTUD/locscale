"""LocScale2 feature-enhance pipeline, self-contained.

Everything runs in memory on numpy arrays; nothing is written to disk except the FDR mask's
temporary files, which the vendored confidence-map code needs. The bundle depends only on
external packages (torch, numpy, scipy, scikit-learn, mrcfile) -- not on the locscale
package.

Stages
------
1.  Mask. If the user supplies one it is used as-is; otherwise an FDR confidence map is
    computed at 1% FDR and returned as an output so it can be inspected and reused.
2.  Preprocess: resample to 1 A/voxel and standardise.
3.  Cube the map, keeping only cubes with signal under the mask.
4.  EMmerNet with Monte-Carlo dropout -> per-cube mean and variance.
5.  Reassemble both, and resample back to the original grid. The mean *is* the
    feature-enhanced map.
6.  Windowed amplitude scaling, GPU-accelerated: target = the input map, reference = the
    feature-enhanced map -> the baseline map.
7.  pVDDT from the three: how far the feature-enhanced map departs from the baseline,
    measured in units of the calibrated standard error.
"""
import os
import tempfile

import numpy as np

from .vendored import mapops
from .emmernet import get_device, load_emmernet, predict_monte_carlo
from .windowed_scaling import run_windowed_scaling


# Palette from LocScale's create_and_save_chimera_script(): blue = the enhanced map sits
# well below the baseline, green = agreement, red = well above.
PVDDT_PALETTE = "-95,#0000ff:-80,#00ffff:0,#00ff00:80,#ffff00:95,#ff0000"


class Cancelled(Exception):
    """Raised by a caller's status/progress callback to abort the run.

    The pipeline never raises or catches this itself; it simply lets it propagate out of
    run_feature_enhance. Cancellation is therefore cooperative and only takes effect where
    the pipeline calls back, i.e. between stages and between EMmerNet batches / scaling
    chunks. The FDR mask stage makes no callbacks and cannot be interrupted.
    """


def _noop(*args, **kwargs):
    pass


def compute_fdr_mask(emmap, apix, window_size=None, fdr=0.01, work_dir=None):
    """FDR confidence map at the given FDR, binarised into a mask.

    Mirrors locscale's run_FDR: window size defaults to 10% of the box, and the raw
    confidence map is thresholded at 0.99.
    """
    if window_size is None:
        window_size = int(round(0.1 * emmap.shape[0]))
        window_size = max(8, window_size)

    work_dir = work_dir or tempfile.mkdtemp(prefix="locscale2_fdr_")
    confidence_map, _ = mapops.compute_FDR_confidenceMap_easy(
        emmap, apix=apix, window_size=window_size, fdr=fdr, folder=work_dir,
        remove_temp_files=True,
    )
    mask = mapops.binarise_map(confidence_map, threshold=0.99,
                               return_type="int", threshold_type="gteq")
    return np.asarray(mask, dtype=np.float32)


def compute_pvddt(feature_enhanced, baseline, variance, n_samples, data_dir=None):
    """Per-voxel confidence, scaled to [-100, +100].

    z = (feature_enhanced - baseline) / calibrated standard error, then the normal CDF
    mapped onto [-100, +100]. The calibrator is the isotonic regression shipped with
    LocScale; NaNs (zero variance) are pushed to +100, i.e. maximally significant, which is
    what the reference does.
    """
    import pickle
    from scipy.stats import norm

    data_dir = data_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    calibrator_path = os.path.join(data_dir, "calibrator_locscale_target_seed_42.pickle")
    with open(calibrator_path, "rb") as handle:
        calibrator = pickle.load(handle)

    standard_error = np.sqrt(np.maximum(variance, 0)) / np.sqrt(n_samples)
    calibrated = calibrator.predict(standard_error.flatten()).reshape(standard_error.shape)

    with np.errstate(divide="ignore", invalid="ignore"):
        z = (feature_enhanced - baseline) / calibrated
    z[~np.isfinite(z)] = 100.0

    return norm.cdf(z) * 200 - 100


def run_feature_enhance(emmap, apix, mask=None, model_type="high_context",
                        monte_carlo_iterations=15, batch_size=8, cube_size=32, stride=16,
                        window_size=25, scaling_chunk=4096, use_gpu=True, gpu_id=None,
                        status_callback=None, progress_callback=None):
    """Run the whole pipeline on in-memory arrays.

    Returns a dict with 'feature_enhanced', 'baseline', 'pvddt', 'variance', and 'mask'
    (the last is None when the caller supplied one). Every array is on the input grid.
    """
    status = status_callback or _noop
    progress = progress_callback or _noop
    emmap = np.asarray(emmap, dtype=np.float32)

    # ---- 1. mask ----------------------------------------------------------
    computed_mask = None
    if mask is None:
        status("Computing FDR confidence mask (no mask supplied)...")
        computed_mask = compute_fdr_mask(emmap, apix)
        mask = computed_mask
    mask = np.asarray(mask, dtype=np.float32)
    if mask.shape != emmap.shape:
        raise ValueError(f"mask shape {mask.shape} does not match map shape {emmap.shape}")

    # ---- 2. preprocess ----------------------------------------------------
    status("Preprocessing (resampling to 1 A/voxel, standardising)...")
    emmap_resampled = mapops.resample_map(emmap, apix=apix, apix_new=1)
    emmap_preprocessed = mapops.standardize_map(emmap_resampled)
    mask_preprocessed = mapops.resample_map(mask, apix=apix, apix_new=1)

    # ---- 3. cube ----------------------------------------------------------
    status("Extracting cubes...")
    cubes_dictionary, cubes_array, signal_cubes = mapops.get_cubes(
        emmap_preprocessed, cube_size=cube_size, step_size=stride, mask=mask_preprocessed)
    status(f"{len(cubes_array)} cubes with signal")

    # ---- 4. EMmerNet, Monte-Carlo dropout ---------------------------------
    device = get_device(prefer_gpu=use_gpu, gpu_id=gpu_id)
    status(f"Loading EMmerNet ({model_type}) on {device}...")
    model, device = load_emmernet(model_type=model_type, device=device)

    status(f"Running {monte_carlo_iterations} Monte-Carlo passes...")
    mean_cubes, var_cubes = predict_monte_carlo(
        cubes_array, model, device, num_samples=monte_carlo_iterations,
        batch_size=batch_size,
        progress_callback=lambda i, n: progress("EMmerNet", i, n))

    # ---- 5. reassemble ----------------------------------------------------
    status("Assembling cubes...")
    feature_enhanced = _assemble(mean_cubes, cubes_dictionary, emmap_preprocessed.shape,
                                 apix, emmap.shape)
    variance = _assemble(var_cubes, cubes_dictionary, emmap_preprocessed.shape,
                         apix, emmap.shape)

    # ---- 6. windowed scaling -> baseline -----------------------------------
    status(f"Windowed amplitude scaling on {device} (window {window_size})...")
    baseline = run_windowed_scaling(
        target_map=emmap, reference_map=feature_enhanced, mask=mask,
        wn=window_size, device=str(device), chunk=scaling_chunk,
        progress_callback=lambda i, n: progress("Scaling", i, n))

    # ---- 7. pVDDT ----------------------------------------------------------
    status("Computing pVDDT...")
    pvddt = compute_pvddt(feature_enhanced, baseline, variance,
                          n_samples=monte_carlo_iterations)

    status("Done.")
    return {
        "feature_enhanced": feature_enhanced,
        "baseline": baseline,
        "pvddt": pvddt,
        "variance": variance,
        "mask": computed_mask,       # None when the caller supplied a mask
    }


def _assemble(cubes, cubes_dictionary, preprocessed_shape, apix, output_shape):
    """Put cubes back in place, then resample to the original grid."""
    filled = mapops.replace_cubes_in_dictionary(cubes, cubes_dictionary)
    assembled = mapops.assemble_cubes(filled, preprocessed_shape, average=True)
    return mapops.resample_map(assembled, apix=1, apix_new=apix, assert_shape=output_shape)
