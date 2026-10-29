"""Two-stage, resumable family calibration and independent evaluation."""
import argparse,datetime,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from biology import generate
from dependence import evaluate,METHODS
ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/'protocol/design.json'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def freeze():
    names=['src/biology.py','src/dependence.py','src/run_dependence.py','tests/test_dependence.py','protocol/design.json']
    hashes={p:sha(ROOT/p) for p in names};target=ROOT/'protocol/dependence_lock.json'
    if target.exists():assert json.loads(target.read_text())['sources']==hashes,'Frozen dependence source changed'
    else:target.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=hashes),indent=2))


def cells(phase):
    d=json.loads(DESIGN.read_text());p=d['dependence'];result=[]
    if phase=='calibration':scenarios=[(s,p['calibration_per_null_cell']) for s in p['null_scenarios']]
    else:scenarios=[(s,p['null_evaluation_per_cell']) for s in p['null_scenarios']]+[(s,p['alternative_evaluation_per_cell']) for s in p['alternative_scenarios']]+[(s,p['stress_evaluation_per_cell']) for s in p['stress_scenarios']]
    for n in d['sample_sizes']:
        for s,(name,count) in enumerate(scenarios):
            for start in range(0,count,p['chunk_size']):result.append((n,s,name,start,min(start+p['chunk_size'],count)))
    return result


def worker(phase,shard,shards):
    freeze();d=json.loads(DESIGN.read_text());p=d['dependence'];folder=ROOT/f'results/dependence_chunks/{phase}';folder.mkdir(parents=True,exist_ok=True)
    for index,(n,s,name,start,end) in enumerate(cells(phase)):
        if index%shards!=shard:continue
        target=folder/f'{name}_{n}_{start:04d}.npz'
        if target.exists():continue
        stats=[];seconds=[];seeds=[];sizes=[];errors=[];gradients=[]
        for rep in range(start,end):
            seed=p[phase+'_seed']+n*1000+s*100000+rep
            try:
                q=evaluate(generate(n,seed,name));stats.append(q['statistics']);seconds.append(q['seconds']);sizes.append(q['pair_sizes']);gradients.append(q['gradient'])
            except Exception as exc:
                stats.append(np.full((len(METHODS),16),np.nan));seconds.append(np.nan);sizes.append(np.zeros(16));gradients.append(np.full(16,np.nan));errors.append(dict(replicate=rep,error=repr(exc)))
            seeds.append(seed)
        np.savez_compressed(target,statistics=stats,seconds=seconds,seeds=seeds,pair_sizes=sizes,gradients=gradients)
        target.with_suffix('.json').write_text(json.dumps(dict(scenario=name,n=n,start=start,end=end,errors=errors,sha256=sha(target)),indent=2))
        print(phase,name,n,start,'failures',len(errors),flush=True)


def calibrate():
    d=json.loads(DESIGN.read_text());p=d['dependence'];thresholds={}
    for n in d['sample_sizes']:
        genes=[];paths=[]
        for scenario in p['null_scenarios']:
            records=[]
            for f in sorted((ROOT/'results/dependence_chunks/calibration').glob(f'{scenario}_{n}_*.npz')):
                meta=json.loads(f.with_suffix('.json').read_text());assert not meta['errors'],meta
                a=np.load(f);records.append(a['statistics'])
            stat=np.concatenate(records);assert len(stat)==p['calibration_per_null_cell'];assert np.isfinite(stat).all()
            q=stat*stat
            genes.append(np.sort(q.max(axis=2),axis=0)[-25])
            paths.append(np.sort(q.reshape(len(q),len(METHODS),4,4).sum(axis=3).max(axis=2),axis=0)[-25])
        thresholds[str(n)]={'gene':np.max(genes,axis=0).tolist(),'module':np.max(paths,axis=0).tolist(),'gene_by_bank':np.asarray(genes).tolist(),'module_by_bank':np.asarray(paths).tolist()}
    target=ROOT/'results/family_thresholds.json'
    payload=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),rank=476,calibration_per_bank=500,averaged_bank_bound=25/501,methods=METHODS,thresholds=thresholds)
    if target.exists():
        old=json.loads(target.read_text());assert old['thresholds']==thresholds
    else:target.write_text(json.dumps(payload,indent=2))
    print('Family thresholds frozen before independent evaluation.',flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--phase',choices=['calibration','evaluation']);ap.add_argument('--shard',type=int,default=0);a=ap.parse_args()
    assert 1<=a.workers<=4;freeze()
    if a.phase:worker(a.phase,a.shard,a.workers);return
    for phase in ['calibration','evaluation']:
        if phase=='evaluation':calibrate()
        jobs=[]
        for shard in range(a.workers):
            log=open(ROOT/f'work/dependence_{phase}_{shard}.log','w',encoding='utf-8')
            proc=subprocess.Popen([sys.executable,__file__,'--phase',phase,'--shard',str(shard),'--workers',str(a.workers)],stdout=log,stderr=subprocess.STDOUT)
            jobs.append((proc,log))
        for proc,log in jobs:
            code=proc.wait();log.close();assert code==0,f'{phase} worker failed'
        print(phase,'completed',flush=True)


if __name__=='__main__':main()
