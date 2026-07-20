#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration test for helical symmetry averaging (symmetrize_map_helical)
against a real, deposited helical cryo-EM reconstruction.

Fixture: EMD-71543, C4 rotational group coincident with the helical axis, 
twist=-15.77 deg, rise=7.87 A, apix=0.83 A/px, 256^3 box.

The map is downloaded on demand the first time this test runs and cached
under tests/test_data/ (gitignored, consistent with the rest of the LocScale
test suite) for subsequent runs.
"""
import os
import unittest

EMDB_ID = "71543"
EMDB_TWIST = -15.77   # degrees, deposited helical twist per asymmetric unit
EMDB_RISE = 7.87      # Angstrom, deposited helical rise per asymmetric unit
EMDB_PG = "C4"        # deposited point group, coincident with the helical axis
EMDB_URL = "https://files.rcsb.org/pub/emdb/structures/EMD-{0}/map/emd_{0}.map.gz".format(EMDB_ID)


class TestHelicalSymmetryRealMap(unittest.TestCase):

    def setUp(self):
        from locscale.utils.file_tools import get_locscale_path
        self.locscale_path = get_locscale_path()
        self.cache_dir = os.path.join(self.locscale_path, "locscale", "tests", "test_data", "helical_cache")
        os.makedirs(self.cache_dir, exist_ok=True)
        self.map_path = os.path.join(self.cache_dir, "emd_{}.map".format(EMDB_ID))
        if not os.path.exists(self.map_path):
            try:
                self._download_map()
            except Exception as exc:
                raise unittest.SkipTest(
                    "Could not download real test map EMD-{} (no network access?): {}".format(EMDB_ID, exc))

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

    def test_helical_symmetrisation_matches_deposited_parameters(self):
        import torch
        from locscale.include.symmetry_emda.symmetrize_map import symmetrize_map_helical
        from locscale.include.emmer.ndimage.map_utils import load_map
        from locscale.include.emmer.ndimage.map_tools import compute_real_space_correlation as rscc

        emmap, apix = load_map(self.map_path)
        print("Loaded EMD-{}: shape {}, apix {}".format(EMDB_ID, emmap.shape, apix))

        # n_steps=3 and complex64 keep this to a couple of minutes on CPU; the
        # automatic (larger) n_steps and complex128 default are for production
        # use, not this correctness check.
        kwargs = dict(n_steps=3, dtype=torch.complex64)

        sym_correct = symmetrize_map_helical(emmap, apix, twist=EMDB_TWIST, rise=EMDB_RISE, pg=EMDB_PG, **kwargs)
        rscc_correct = rscc(emmap, sym_correct)
        print("RSCC with deposited helical parameters (twist={}, rise={}, pg={}): {:.4f}".format(
            EMDB_TWIST, EMDB_RISE, EMDB_PG, rscc_correct))

        # A genuine helical reconstruction already closely obeys its own
        # symmetry, so averaging over the correct operators should reinforce
        # signal rather than wash it out.
        self.assertGreater(rscc_correct, 0.9)

        # Negative control: the wrong handedness (twist sign flipped) is a
        # different, incompatible symmetry -- correlation should drop sharply.
        sym_wrong_hand = symmetrize_map_helical(emmap, apix, twist=-EMDB_TWIST, rise=EMDB_RISE, pg=EMDB_PG, **kwargs)
        rscc_wrong_hand = rscc(emmap, sym_wrong_hand)
        print("RSCC with flipped twist handedness: {:.4f}".format(rscc_wrong_hand))
        self.assertLess(rscc_wrong_hand, rscc_correct - 0.25)

        # Negative control: the wrong point-group order (C3 instead of the
        # deposited C4) should correlate distinctly worse. 
        sym_wrong_pg = symmetrize_map_helical(emmap, apix, twist=EMDB_TWIST, rise=EMDB_RISE, pg="C3", **kwargs)
        rscc_wrong_pg = rscc(emmap, sym_wrong_pg)
        print("RSCC with pg=C3 instead of deposited C4: {:.4f}".format(rscc_wrong_pg))
        self.assertLess(rscc_wrong_pg, rscc_correct - 0.2)


if __name__ == "__main__":
    unittest.main()
