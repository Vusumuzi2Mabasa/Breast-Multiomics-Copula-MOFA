import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from simulation import generate
from eqtl import bh,nested,analyse
from integration import prep_fit,prep_apply,mofa,project


class Science(unittest.TestCase):
    def test_replay_and_layers(self):
        a=generate(120,81,'subtype_eqtl');b=generate(120,81,'subtype_eqtl')
        np.testing.assert_equal(a['counts'],b['counts']);self.assertNotIn('protein',a)
        self.assertTrue(np.isin(a['genotype'],[0,1,2]).all());self.assertTrue(np.isfinite(a['RNA']).all())
        self.assertEqual(a['genotype'].shape,(120,24));self.assertEqual(a['RNA'].shape,(120,12))

    def test_bh_known_values(self):
        np.testing.assert_allclose(bh([.01,.04,.03,.2]),[.04,.053333333333,.053333333333,.2])

    def test_nested_invariance(self):
        r=np.random.default_rng(1);x=r.normal(size=(120,3));y=x[:,2]+r.normal(size=120)
        a=np.c_[np.ones(120),x[:,:2]];b=np.c_[a,x[:,2]]
        p,stat,_=nested(a,b,y);pp,ss,_=nested(a,b,3+5*y)
        self.assertAlmostEqual(stat,ss,10);self.assertAlmostEqual(p,pp,12)

    def test_truth_and_sign(self):
        d=generate(1500,187,'subtype_eqtl');rows=analyse(d)
        self.assertEqual(sum(r['true_eqtl'] for r in rows),9);self.assertEqual(sum(r['true_interaction'] for r in rows),6)
        self.assertGreater(rows[0]['causal_interaction_estimate'],0);self.assertLess(rows[3]['causal_interaction_estimate'],0)
        self.assertTrue(all(not r['true_eqtl'] for r in analyse(generate(100,54,'abundance_only'))))

    def test_training_preparation(self):
        d=generate(100,43,'shared_eqtl');views=[d['genotype'],d['RNA']];a,p=prep_fit(views,d['covariates']);b=prep_apply(views,d['covariates'],p)
        for x,y in zip(a,b):np.testing.assert_allclose(x,y)
        # Applying a shift to test data does not alter training transformation parameters.
        before=[x.copy() for x,_ in p];prep_apply([v+50 for v in views],d['covariates'],p)
        for x,(y,_) in zip(before,p):np.testing.assert_equal(x,y)

    def test_official_mofa_and_projection(self):
        d=generate(100,166,'subtype_eqtl');v,_=prep_fit([d['genotype'],d['RNA']],d['covariates'])
        a=mofa(v,2,71,35);b=mofa(v,2,72,35,True)
        self.assertNotEqual(a['initialization_hash'],b['initialization_hash'])
        z=project(a,v);self.assertEqual(z.shape,(100,2));self.assertTrue(np.isfinite(z).all())


if __name__=='__main__':unittest.main()
