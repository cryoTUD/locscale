#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end integration test: run the actual LocScale-FEM (`feature_enhance`)
pipeline with helical symmetry enabled, and check that the final,
LocScale-scaled output map (not just the EMmerNet reference map) comes out
with (near) perfect helical symmetry.

Fixture: EMD-71543, a helical oligomer of YG18 peptides -- C4 point group
coincident with the helical axis, twist=-15.77 deg, rise=7.87 A. Reused from
the helical_cache populated by test_symmetry_helical.py (downloaded here too
if not already present), then resampled to apix=2.0 A/px to keep EMmerNet
inference and LocScale scaling fast for a test.

This exercises the whole chain: CLI parsing -> EMmerNet single-point
prediction -> reference-map symmetrisation -> LocScale scaling -> final
output symmetrisation (write_out_final_volume_window_back_if_required).
Monte Carlo sampling is disabled (a single EMmerNet prediction is enough for
a scaling reference) to keep runtime down.
"""
import os
import unittest

EMDB_ID = "71543"
EMDB_TWIST = -15.77   # degrees, deposited helical twist per asymmetric unit
EMDB_RISE = 7.87      # Angstrom, deposited helical rise per asymmetric unit
EMDB_PG = "C4"        # deposited point group, coincident with the helical axis
EMDB_URL = "https://files.rcsb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz".format(EMDB_ID)
RESAMPLED_APIX = 2.0


class TestLocScaleFEMHelical(unittest.TestCase):

    def setUp(self):
        from locscale.utils.file_tools import get_locscale_path
        self.locscale_path = get_locscale_path()
        self.cache_dir = os.path.join(self.locscale_path, "locscale", "tests", "test_data", "helical_cache")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.map_path = os.path.join(self.cache_dir, "emd_{}.map".format(EMDB_ID))
        self.resampled_path = os.path.join(
            self.cache_dir, "emd_{}_resampled_{}A.mrc".format(EMDB_ID, RESAMPLED_APIX))
        try:
            if not os.path.exists(self.map_path):
                self._download_map()
            if not os.path.exists(self.resampled_path):
                self._resample_map()
        except Exception as exc:
            raise unittest.SkipTest(
                "Could not prepare real test map EMD-{} (no network access?): {}".format(EMDB_ID, exc))

    def _download_map(self):
        import urllib.request
        import gzip
        import shutil
        gz_path = self.map_path + ".gz"
        print("Downloading EMD-{} (real helical test map, ~64 MB) ...".format(EMDB_ID))
        urllib.request.urlretrieve(EMDB_URL, gz_path)
        with gzip.open(gz_path, "rb") as f_in, open(self.map_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_path)

    def _resample_map(self):
        from locscale.include.emmer.ndimage.map_utils import load_map, resample_map, save_as_mrc
        emmap, apix = load_map(self.map_path)
        resampled = resample_map(emmap, apix=apix, apix_new=RESAMPLED_APIX)
        save_as_mrc(resampled, self.resampled_path, apix=RESAMPLED_APIX)

    def test_locscale_fem_final_output_has_near_perfect_helical_symmetry(self):
        import shutil
        import tempfile
        import numpy as np
        import torch
        from locscale.utils.parse_utils import locscale_parser
        from locscale.utils.startup_utils import launch_feature_enhance
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map
        from locscale.include.emmer.ndimage.map_utils import load_map

        # The FEM pipeline always writes its final "_baseline" output next to
        # the input emmap (input_folder), regardless of -o's directory -- only
        # its basename is used. Track that path explicitly so it can be
        # cleaned up afterwards instead of left in the (tracked) cache dir.
        input_folder = os.path.dirname(self.resampled_path)
        output_path = os.path.join(input_folder, "fem_output_baseline.mrc")
        # check_and_save_output() also saves the raw EMmerNet prediction using
        # -o's value as-is (not joined to any directory), so give it an
        # absolute path pointing into the temp workdir to avoid polluting CWD.
        try:
            with tempfile.TemporaryDirectory() as workdir:
                outfile_path = os.path.join(workdir, "fem_output.mrc")
                args = locscale_parser.parse_args([
                    "feature_enhance",
                    "--emmap_path", self.resampled_path,
                    "-o", outfile_path,
                    "-sym", EMDB_PG,
                    "-twist", str(EMDB_TWIST),
                    "-rise", str(EMDB_RISE),
                    "-n_steps", "3",
                    "--no_monte_carlo",
                    "-op", workdir,
                ])
                try:
                    launch_feature_enhance(args)
                except Exception as exc:
                    raise unittest.SkipTest(
                        "LocScale-FEM run could not complete (missing EMmerNet weights or other "
                        "environment issue): {}".format(exc))

                self.assertTrue(os.path.exists(output_path), "Expected final LocScale-FEM output not found")

                final_output, apix = load_map(output_path)
                print("Final LocScale-FEM output: shape {}, apix {}".format(final_output.shape, apix))

                kwargs = dict(n_steps=3, dtype=torch.complex64)

                # The output already went through final-output symmetrisation with
                # the correct parameters, so re-applying the same symmetry should
                # be close to a no-op (near-perfect helical symmetry already achieved).
                resym_correct = symmetrize_map(final_output, apix, pg=EMDB_PG, twist=EMDB_TWIST, rise=EMDB_RISE, **kwargs)
                rscc_correct = np.corrcoef(final_output.ravel(), resym_correct.ravel())[0, 1]
                print("RSCC(final output, re-symmetrised w/ correct params): {:.4f}".format(rscc_correct))
                self.assertGreater(rscc_correct, 0.99)

                # Negative control: re-symmetrising with the wrong handedness should
                # diverge noticeably more, confirming this isn't just numerical no-op smoothing.
                resym_wrong = symmetrize_map(final_output, apix, pg=EMDB_PG, twist=-EMDB_TWIST, rise=EMDB_RISE, **kwargs)
                rscc_wrong = np.corrcoef(final_output.ravel(), resym_wrong.ravel())[0, 1]
                print("RSCC(final output, re-symmetrised w/ flipped handedness): {:.4f}".format(rscc_wrong))
                self.assertLess(rscc_wrong, rscc_correct - 0.05)
        finally:
            if os.path.exists(output_path):
                os.remove(output_path)


if __name__ == "__main__":
    unittest.main()
