import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from metrics import full_subspace_recovery


class TestRecovery(unittest.TestCase):
    def test_missing_dimension_is_zero(self):
        rng=np.random.default_rng(12);truth=rng.normal(size=(100,3))
        self.assertAlmostEqual(full_subspace_recovery(truth[:,:2],truth),2/3,places=12)

    def test_rotation_and_sign_invariance(self):
        rng=np.random.default_rng(13);truth=rng.normal(size=(100,3))
        rotation=np.linalg.qr(rng.normal(size=(3,3)))[0]
        self.assertAlmostEqual(full_subspace_recovery(truth@rotation,truth),1,places=12)

    def test_no_recovered_dimension(self):
        rng=np.random.default_rng(14);truth=rng.normal(size=(100,3))
        self.assertEqual(full_subspace_recovery(np.ones((100,3)),truth),0)


if __name__=='__main__':unittest.main()
