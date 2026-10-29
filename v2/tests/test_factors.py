import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biology import generate
from factors import fit_preprocess,apply_preprocess,gibbs,official_mofa
from integration import datasets


class TestIntegration(unittest.TestCase):
    def test_distinct_official_initializations(self):
        d=generate(100,713,'both');v,_=fit_preprocess(d['views'],d['covariates'])
        a=official_mofa(v,2,811,iterations=30,perturb=False)
        b=official_mofa(v,2,812,iterations=30,perturb=True)
        self.assertNotEqual(a['initialization_hash'],b['initialization_hash'])
        self.assertEqual(b['initialization'],'PCA plus jitter')

    def test_independent_country_split(self):
        p=dict(test_size=60,unseen_test_size=40)
        a,b,c=datasets(100,72105,'both',p)
        self.assertNotIn(4,a['country']);self.assertNotIn(4,b['country']);self.assertTrue(np.all(c['country']==4))
        self.assertNotEqual(a['seed'],b['seed']);self.assertNotEqual(a['seed'],c['seed'])

    def test_training_preprocessing_and_binary_view(self):
        d=generate(120,711,'missing_coordination');v,params=fit_preprocess(d['views'],d['covariates'])
        np.testing.assert_equal(v[1],d['somatic'])
        rebuilt=apply_preprocess(d['views'],d['covariates'],params)
        for a,b in zip(v,rebuilt):np.testing.assert_allclose(a,b,equal_nan=True)

    def test_rank_sampler_monotone_invariance(self):
        rng=np.random.default_rng(712);x=rng.normal(size=(80,6));x[:,0]=(x[:,0]>0).astype(float)
        config=dict(gibbs_burn=50,gibbs_draws=100,gibbs_thin=2)
        a=gibbs(x,2,931,True,config);b=gibbs(np.exp(x),2,931,True,config)
        np.testing.assert_allclose(a['correlation'],b['correlation'],atol=1e-12)
        self.assertGreater(np.linalg.eigvalsh(a['correlation']).min(),0)


if __name__=='__main__':unittest.main()
