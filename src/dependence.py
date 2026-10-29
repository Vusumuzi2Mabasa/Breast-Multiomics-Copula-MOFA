"""Independent calibration of a prespecified RNA/protein dependence test.

No feature search: one known pair; no use of factor-model fits or TNBC means.
This companion experiment tests a different estimand from factor integration.
"""
from pathlib import Path
import json
import numpy as np
from scipy.special import ndtr,ndtri
from scipy.stats import nbinom,rankdata
ROOT=Path(__file__).resolve().parents[1]


def statistic(n,seed,rho):
    rng=np.random.default_rng(seed);y=rng.binomial(1,.3,n)
    a=rng.normal(size=n);b=rho*y*a+np.sqrt(1-(rho*y)**2)*rng.normal(size=n)
    rna=np.log2(nbinom.ppf(np.clip(ndtr(a),1e-8,1-1e-8),8,8/(8+80))+1)
    protein=2+.7*b
    normal=lambda v:ndtri((rankdata(v)-.5)/len(v))
    vals=[]
    for x,z in [(rna,protein),(normal(rna),normal(protein))]:
        r=[np.corrcoef(x[y==g],z[y==g])[0,1] for g in [0,1]]
        vals.append(abs(np.arctanh(np.clip(r[1],-.999,.999))-np.arctanh(np.clip(r[0],-.999,.999)))/np.sqrt(1/(sum(y==0)-3)+1/(sum(y==1)-3)))
    return vals


def main():
    out=ROOT/'results/dependence';out.mkdir(parents=True,exist_ok=True)
    records=[];thresholds={}
    for n in [180,360]:
        cal=np.array([statistic(n,39102901+n*1000+r,0) for r in range(500)])
        cut=np.sort(cal,axis=0)[-25]
        thresholds[str(n)]=dict(thresholds=cut.tolist(),calibration_n=500,
                                strict_rank=476,averaged_null_bound=25/501)
        np.savez_compressed(out/f'calibration_{n}.npz',statistics=cal)
        for name,rho in [('null',0),('dependence_change',.65)]:
            bank=0 if name=='null' else 1000000
            for r in range(500):
                vals=statistic(n,49102901+bank+n*1000+r,rho)
                for m,(score,threshold) in enumerate(zip(vals,cut)):
                    records.append(dict(n=n,scenario=name,replicate=r,method=['Pearson log','Copula normal scores'][m],statistic=float(score),threshold=float(threshold),reject=bool(score>threshold)))
    (out/'thresholds.json').write_text(json.dumps(thresholds,indent=2))
    (out/'records.json').write_text(json.dumps(records))
    print('Finished 1000 calibration and 2000 evaluation datasets; 4000 method records.')


if __name__=='__main__':main()
