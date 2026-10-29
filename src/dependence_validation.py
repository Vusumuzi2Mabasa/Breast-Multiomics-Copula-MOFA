"""Prospective, uniform calibration extension; preserve all initial outputs."""
from pathlib import Path
import hashlib,json,datetime
import numpy as np
from dependence import statistic
ROOT=Path(__file__).resolve().parents[1]


def main():
    protocol=json.loads((ROOT/'protocol/dependence_validation.json').read_text())
    lockfile=ROOT/'protocol/dependence_validation_lock.json'
    names=['protocol/dependence_validation.json','src/dependence_validation.py','src/dependence.py']
    hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    if lockfile.exists():assert json.loads(lockfile.read_text())['sources']==hashes
    else:lockfile.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=hashes),indent=2))
    out=ROOT/'results/dependence_validation';out.mkdir(exist_ok=True)
    records=[];thresholds={}
    for n in protocol['sample_sizes']:
        cal=np.array([statistic(n,protocol['calibration_seed']+n*1000+r,0) for r in range(5000)])
        cut=np.sort(cal,axis=0)[-250]
        thresholds[str(n)]=dict(thresholds=cut.tolist(),calibration_n=5000,strict_rank=4751,averaged_null_bound=250/5001)
        np.savez_compressed(out/f'calibration_{n}.npz',statistics=cal)
        for cell,rho in enumerate(protocol['rho_alternatives']):
            name='null' if rho==0 else f'rho_{rho:.2f}'
            for r in range(1000):
                vals=statistic(n,protocol['evaluation_seed']+cell*1000000+n*1000+r,rho)
                for method,(score,threshold) in enumerate(zip(vals,cut)):
                    records.append(dict(n=n,scenario=name,rho=rho,replicate=r,method=['Pearson log','Copula normal scores'][method],statistic=float(score),threshold=float(threshold),reject=bool(score>threshold)))
            print(n,name,'complete',flush=True)
    (out/'thresholds.json').write_text(json.dumps(thresholds,indent=2))
    (out/'records.json').write_text(json.dumps(records))
    print('Finished 10000 calibration and 8000 fresh evaluation datasets.')


if __name__=='__main__':main()
