import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from biology import generate
from dependence import evaluate


class TestBiologicalControls(unittest.TestCase):
    def test_genotype_coding_invariance(self):
        d=generate(240,761,'both');a=evaluate(d)['statistics']
        e=dict(d);e['germline']=2-d['germline'];b=evaluate(e)['statistics']
        np.testing.assert_allclose(a,b,atol=1e-6,rtol=1e-6)

    def test_additive_subtype_shift_removed_by_conditional_margins(self):
        d=generate(240,762,'null');a=evaluate(d)['statistics']
        e=dict(d);e['RNA']=d['RNA']+2*d['y'][:,None];e['protein']=d['protein']-.5*d['y'][:,None]
        b=evaluate(e)['statistics']
        # Raw-rank partial Spearman is not affine-invariant across subgroups;
        # conditional location/scale statistics should remove this exact shift.
        np.testing.assert_allclose(a[2:],b[2:],atol=1e-6,rtol=1e-6)


if __name__=='__main__':unittest.main()
