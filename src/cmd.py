"""
Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the
program 'LocScale2' written by the Author(s).

Copyright (C) 2026 Alok Bharadwaj and Arjen J. Jakobi

The `locscale2` ChimeraX command.
"""
import numpy as np
from chimerax.core.commands import BoolArg, CmdDesc, EnumOf, IntArg
from chimerax.map import MapArg, volume_from_grid_data
from chimerax.map_data import ArrayGridData

from .emmernet import available_models
from .pipeline import PVDDT_PALETTE, run_feature_enhance


def show_volume(session, array, template, name, show=True):
    """Open a numpy array as a volume on the template's grid, so it stays superimposed."""
    grid = ArrayGridData(np.ascontiguousarray(array, dtype=np.float32),
                         origin=template.data.origin, step=template.data.step, name=name)
    grid.name = name
    volume = volume_from_grid_data(grid, session)
    volume.display = show
    if show:
        # Setting display only *schedules* the surface to be built on the next frame, so a
        # command issued straight afterwards -- `color sample`, below -- would find no
        # surface on this model and fail. Build it now. (ChimeraX does this itself only
        # when session.in_script is set, which is not the case for a tool or a command.)
        volume.update_drawings()
    return volume


def show_results(session, results, template):
    """Open the outputs and colour the feature-enhanced map by pVDDT."""
    opened = {}

    if results.get("mask") is not None:
        # only present when we computed it; shown off so it does not obscure the maps
        opened["mask"] = show_volume(session, results["mask"], template,
                                     "LocScale2 FDR mask", show=False)

    opened["feature_enhanced"] = show_volume(session, results["feature_enhanced"], template,
                                             "LocScale2 feature enhanced")
    opened["baseline"] = show_volume(session, results["baseline"], template,
                                     "LocScale2 baseline", show=False)
    # pVDDT is a colour source rather than something to look at directly
    opened["pvddt"] = show_volume(session, results["pvddt"], template,
                                  "LocScale2 pVDDT", show=False)

    from chimerax.core.commands import run
    run(session, "color sample #{} map #{} palette {}".format(
        opened["feature_enhanced"].id_string, opened["pvddt"].id_string, PVDDT_PALETTE))
    session.logger.info(
        "LocScale2: coloured the feature-enhanced map by pVDDT "
        "(blue = below baseline, green = agreement, red = above).")

    return opened


def show_locscale_result(session, results, template):
    """Open the amplitude-scaled map (and the FDR mask, if one was computed)."""
    opened = {}
    if results.get("mask") is not None:
        opened["mask"] = show_volume(session, results["mask"], template,
                                     "LocScale2 FDR mask", show=False)
    opened["locscale"] = show_volume(session, results["locscale"], template,
                                     "LocScale amplitude scaled")
    session.logger.info("LocScale2: amplitude scaling done (no feature enhancement).")
    return opened


def locscale2(session, inputMap=None, inputMask=None, model="high_context",
              monteCarloIterations=15, batchSize=8, windowSize=25, cubeSize=32, stride=16,
              scalingChunk=4096, useGpu=True, gpuId=None):
    """Feature-enhance a cryo-EM map with EMmerNet and colour it by pVDDT.

    Runs synchronously and blocks the UI; the GUI tool runs the same pipeline on a thread.
    """
    if inputMap is None:
        raise ValueError("Provide an input map, e.g.: locscale2 #1")

    emmap = inputMap.data.full_matrix()
    step = inputMap.data.step
    apix = float(step[0])
    if not all(abs(s - apix) < 1e-4 for s in step):
        session.logger.warning(
            "Non-isotropic voxel size {}; using {:.4f} A for all axes.".format(step, apix))

    mask = None
    if inputMask is not None:
        mask = inputMask.data.full_matrix()
        if mask.shape != emmap.shape:
            raise ValueError("Mask shape {} does not match map shape {}.".format(
                mask.shape, emmap.shape))
    else:
        session.logger.info("LocScale2: no mask given, an FDR mask will be computed.")

    results = run_feature_enhance(
        emmap=emmap, apix=apix, mask=mask, model_type=model,
        monte_carlo_iterations=monteCarloIterations, batch_size=batchSize,
        cube_size=cubeSize, stride=stride, window_size=windowSize,
        scaling_chunk=scalingChunk, use_gpu=useGpu, gpu_id=gpuId,
        status_callback=lambda m: session.logger.info("LocScale2: " + m),
    )
    return show_results(session, results, inputMap)


locscale2_desc = CmdDesc(
    # MapArg, not the Volume class: `usage` needs a real Annotation, and a bare model class
    # makes it fail with "type object 'Volume' has no attribute 'url'".
    required=[("inputMap", MapArg)],
    keyword=[
        ("inputMask", MapArg),
        ("model", EnumOf(available_models())),
        ("monteCarloIterations", IntArg),
        ("batchSize", IntArg),
        ("windowSize", IntArg),
        ("cubeSize", IntArg),
        ("stride", IntArg),
        ("scalingChunk", IntArg),
        ("useGpu", BoolArg),
        ("gpuId", IntArg),
    ],
    synopsis="Feature-enhance a cryo-EM map with EMmerNet and colour it by pVDDT",
)
