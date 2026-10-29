import sys, unittest
from pathlib import Path
import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from models import rank_metadata,pack_metadata,latent_update_packed,assert_order,fit_copula,normal_scores_fit,normal_scores_apply,split_rhat
from simulate import generate


class TestModels(unittest.TestCase):
    def test_ordering_ties_and_missing(self):
        rng=np.random.default_rng(81);x=rng.normal(size=(100,8))
        x[:,0]=rng.binomial(1,.1,100);x[:,1]=rng.binomial(2,.2,100)
        x[::4,2]=np.nan
        meta=rank_metadata(x);packed=pack_metadata(meta,*x.shape)
        z=np.zeros_like(x)
        for j,(order,_,_) in enumerate(meta):z[order,j]=ndtri((rankdata(x[order,j])-.5)/len(order))
        for _ in range(100):
            z=latent_update_packed(z,np.zeros_like(z),packed,rng)
            assert_order(z,meta)
        self.assertTrue(np.isfinite(z).all())
        self.assertGreater(np.std(z[np.isnan(x)]),.5)

    def test_monotone_invariance(self):
        rng=np.random.default_rng(98);x=rng.normal(size=(100,4))
        a,m=normal_scores_fit(x);b,_=normal_scores_fit(np.exp(x))
        np.testing.assert_allclose(a,b)
        np.testing.assert_allclose(normal_scores_apply(x,m),a)

    def test_reproducible_generator_and_panels(self):
        a=generate(200,151,'missing');b=generate(200,151,'missing')
        for x,y in zip(a['views'],b['views']):np.testing.assert_equal(x,y)
        self.assertEqual([v.shape[1] for v in a['views']],[6,6,24,12])
        self.assertTrue(set(np.unique(a['views'][0])).issubset({0,1,2}))
        self.assertTrue(set(np.unique(a['views'][1])).issubset({0,1}))
        self.assertGreater(np.isnan(a['views'][3]).mean(),.05)

    def test_known_correlation(self):
        rng=np.random.default_rng(775);x=rng.multivariate_normal([0,0],[[1,.6],[.6,1]],size=600)
        x[:,0]=np.exp(x[:,0])
        result=fit_copula(x,1,71,burn=300,draws=600)
        self.assertLess(abs(result['correlation'][0,1]-.6),.12)
        np.testing.assert_allclose(np.diag(result['correlation']),1)
        self.assertGreater(np.linalg.eigvalsh(result['correlation']).min(),0)

    def test_rhat_detects_shift(self):
        rng=np.random.default_rng(19);x=rng.normal(size=(2,500,4))
        self.assertLess(split_rhat(x).max(),1.03)
        x[1]+=3
        self.assertGreater(split_rhat(x).min(),1.4)


if __name__=='__main__':unittest.main()
