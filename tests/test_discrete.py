"""Check binary latent correlation against an independent likelihood fit."""
import sys,unittest
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import ndtr,ndtri
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from models import fit_copula


class TestBinaryCopula(unittest.TestCase):
    def test_against_tetrachoric_likelihood(self):
        rng=np.random.default_rng(15521)
        z=rng.multivariate_normal([0,0],[[1,.6],[.6,1]],size=1000)
        x=(z>ndtri([.7,.8])).astype(float)
        counts=np.array([sum((x[:,0]==a)&(x[:,1]==b)) for a,b in [(0,0),(0,1),(1,0),(1,1)]])
        pa=np.mean(x[:,0]==0);pb=np.mean(x[:,1]==0)
        a=ndtri(pa);b=ndtri(pb)
        def objective(rho):
            p00=quad(lambda v:np.exp(-v*v/2)/np.sqrt(2*np.pi)*ndtr((b-rho*v)/np.sqrt(1-rho*rho)),-np.inf,a,epsabs=1e-10)[0]
            probs=np.array([p00,pa-p00,pb-p00,1-pa-pb+p00])
            return -np.sum(counts*np.log(np.maximum(probs,1e-12)))
        mle=minimize_scalar(objective,bounds=(-.9,.9),method='bounded').x
        fit=fit_copula(x,1,812,burn=600,draws=1200)
        self.assertLess(abs(fit['correlation'][0,1]-mle),.12)
        self.assertLess(abs(mle-.6),.12)


if __name__=='__main__':unittest.main()
