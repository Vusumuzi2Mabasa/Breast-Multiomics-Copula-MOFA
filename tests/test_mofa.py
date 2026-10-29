"""Integration check against the actual MOFA+ training engine."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from models import fit_mofa
from simulate import generate,standardise
from benchmark import project_mofa


class TestOfficialMOFA(unittest.TestCase):
    def test_engine_and_new_participant_projection(self):
        train=generate(120,78155,'missing');test=generate(30,78156,'missing')
        a,b=standardise(train['views'],test['views'])
        result=fit_mofa(a,3,491,iterations=150)
        scores=project_mofa(result,b)
        self.assertEqual(scores.shape,(30,3))
        self.assertTrue(np.isfinite(scores).all())
        self.assertGreater(len(result['elbo']),2)
        self.assertTrue(np.isfinite(result['elbo']).all())


if __name__=='__main__':unittest.main()
