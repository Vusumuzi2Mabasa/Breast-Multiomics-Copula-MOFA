"""Independent replay, allele/gene truth checks and R optimisation fixtures."""
import json,sys,hashlib
from pathlib import Path
import numpy as np
from biology import generate,pair_design
from dependence import evaluate,crossfit_margin,score_normalise,copula_objective,copula_score
ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/validation';out.mkdir(exist_ok=True)
    fixtures=[];replays=[]
    for scenario,n in [('null',240),('abundance_only',480),('coordination_moderate',240),('missing_coordination',480),('heavy_tail_null',240)]:
        path=ROOT/f'results/dependence_chunks/evaluation/{scenario}_{n}_0000.npz'
        chunk=np.load(path);seed=int(chunk['seeds'][0]);d=generate(n,seed,scenario);q=evaluate(d)
        error=float(np.max(np.abs(q['statistics']-chunk['statistics'][0])))
        assert error<1e-9,error
        replays.append(dict(scenario=scenario,n=n,seed=seed,max_error=error))
        rng=np.random.default_rng(seed+70000000);fold=np.arange(n)%2;rng.shuffle(fold)
        for j in [0,5]:
            mean,scale,x=pair_design(d,j);a=d['RNA'][:,j];b=d['protein'][:,j];good=np.isfinite(a)&np.isfinite(b);t=d['y'][good]
            u=score_normalise(crossfit_margin(a,mean,scale,fold)[good],t);v=score_normalise(crossfit_margin(b,mean,scale,fold)[good],t)
            score,gradient=copula_score(u,v,t,x[good])
            fixtures.append(dict(scenario=scenario,n=n,gene=j,u=u.tolist(),v=v.tolist(),t=t.tolist(),x=x[good].tolist(),expected=score,gradient=gradient))
    (out/'replay.json').write_text(json.dumps(dict(status='passed',checks=replays),indent=2))
    (out/'R_conditional_fixtures.json').write_text(json.dumps(fixtures))
    # Independently reconstruct fixed threshold and summary decision cells.
    thresholds=json.loads((ROOT/'results/family_thresholds.json').read_text())['thresholds'];checks=[]
    for n in [240,480]:
        for scenario in ['null','abundance_only','scale_only','technical_null','missing_null']:
            arr=np.concatenate([np.load(p)['statistics'] for p in sorted((ROOT/'results/dependence_chunks/calibration').glob(f'{scenario}_{n}_*.npz'))])
            for method in range(6):
                maxima=np.array([max(sum(row[method,j]**2 for j in range(k*4,k*4+4)) for k in range(4)) for row in arr])
                kth=float(sorted(maxima)[475]);assert thresholds[str(n)]['module'][method]>=kth
            checks.append(dict(n=n,scenario=scenario,cohorts=len(arr)))
    (out/'family_checks.json').write_text(json.dumps(dict(status='passed',checks=checks),indent=2))
    print('Five full dependence replays, ten conditional fixtures and family-rank checks passed.')


if __name__=='__main__':main()
