import unittest,sys
from pathlib import Path
import numpy as np
from scipy.optimize._numdiff import approx_derivative
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biology import generate,pair_design
from dependence import copula_objective,copula_score,evaluate


class TestConditionalDependence(unittest.TestCase):
    def test_analytic_gradient(self):
        rng=np.random.default_rng(541);x=np.c_[np.ones(200),rng.normal(size=(200,3))]
        u=rng.normal(size=200);v=.3*u+rng.normal(size=200);theta=np.array([.2,.1,-.1,.05])
        fun=lambda p:copula_objective(p,u,v,x)[0]
        numeric=approx_derivative(fun,theta,method='3-point').ravel()
        np.testing.assert_allclose(numeric,copula_objective(theta,u,v,x)[1],rtol=1e-6,atol=1e-7)

    def test_symmetry_and_sign(self):
        rng=np.random.default_rng(542);x=np.c_[np.ones(300),rng.normal(size=(300,2))];t=rng.binomial(1,.3,300)
        u=rng.normal(size=300);v=.5*t*u+np.sqrt(1-(.5*t)**2)*rng.normal(size=300)
        a,_=copula_score(u,v,t,x);b,_=copula_score(v,u,t,x);c,_=copula_score(u,-v,t,x)
        self.assertAlmostEqual(a,b,places=8);self.assertAlmostEqual(a,-c,places=8)

    def test_gene_modules_and_matching(self):
        d=generate(240,543,'both');e=generate(240,543,'both')
        for a,b in zip(d['views'],e['views']):np.testing.assert_equal(a,b)
        self.assertEqual([a.shape[1] for a in d['views']],[8,8,16,16])
        self.assertEqual(np.flatnonzero(d['beta']).tolist(),[0,1,2,3])
        self.assertEqual(d['counts'].shape,(240,16))
        self.assertTrue(np.isfinite(d['RNA']).all())

    def test_abundance_is_not_coordination_truth(self):
        d=generate(240,544,'abundance_only');self.assertTrue(np.all(d['beta']==0))
        q=evaluate(d);self.assertEqual(q['statistics'].shape,(6,16))
        self.assertTrue(np.isfinite(q['statistics']).all())


if __name__=='__main__':unittest.main()
