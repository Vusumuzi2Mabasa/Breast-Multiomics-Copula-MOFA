"""Resumable, label-free model selection and independent predictive assessment."""
import argparse,datetime,hashlib,json,subprocess,sys
from pathlib import Path
import numpy as np
from integration import datasets,run
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def freeze():
    files=['v2/src/biology.py','v2/src/factors.py','v2/src/integration.py','v2/src/run_integration.py','v2/protocol/design.json','v2/protocol/integration_amendment.json','v2/protocol/initialization_amendment.json','src/models.py','src/benchmark.py','src/metrics.py']
    obj={p:digest(BASE/p) for p in files};target=ROOT/'protocol/integration_lock_r2.json'
    if target.exists():assert json.loads(target.read_text())['sources']==obj,'Frozen integration source changed'
    else:target.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=obj),indent=2))


def export_r(stem,data,folder):
    folder.mkdir(parents=True,exist_ok=True)
    for label,d in zip(['train','test','unseen'],data):
        for name,a in zip(['germline','somatic','RNA','protein'],d['views']):np.savetxt(folder/f'{label}_{name}.csv',a,delimiter=',')
        np.savetxt(folder/f'{label}_covariates.csv',d['covariates'],delimiter=',')
        np.savetxt(folder/f'{label}_y.csv',d['y'],delimiter=',',fmt='%d')


def worker(shard,shards,pilot=False):
    d=json.loads((ROOT/'protocol/design.json').read_text());p=d['integration']
    p.update({k:v for k,v in json.loads((ROOT/'protocol/integration_amendment.json').read_text()).items() if k=='factor_grid'})
    if not pilot:freeze()
    folder=ROOT/('work/integration_pilot' if pilot else 'results/integration');folder.mkdir(exist_ok=True)
    scenarios=['both'] if pilot else p['scenarios'];ns=[240] if pilot else d['sample_sizes'];count=1 if pilot else p['replicates']
    for s,scenario in enumerate(scenarios):
        for n in ns:
            for rep in range(count):
                if rep%shards!=shard:continue
                stem=f'{scenario}_{n}_{rep:02d}';target=folder/(stem+'.json')
                if target.exists():continue
                seed=(101030000 if pilot else p['seed'])+s*100000+n*1000+rep
                data=datasets(n,seed,scenario,p);q=run(*data,seed,p)
                payload=dict(scenario=scenario,n=n,replicate=rep,seed=seed,records=q['records'],predictions=q['predictions'],tunings=q['tunings'],test_y=data[1]['y'].tolist(),unseen_y=data[2]['y'].tolist())
                target.write_text(json.dumps(payload,indent=2,allow_nan=False))
                if rep<2:
                    fixture=q['fixtures'];fixture.update(train_y=data[0]['y'],train_h=data[0]['h'],train_country=data[0]['country'],train_data=np.concatenate(data[0]['views'],axis=1))
                    np.savez_compressed(folder/(stem+'.npz'),**fixture)
                # Raw independent data exported only to ignored work space for
                # the separate official R DIABLO comparator.
                export_r(stem,data,ROOT/'work/diablo_inputs'/stem)
                print(stem,[(r['method'],r['status'],round(r.get('auc',0),3),r['k']) for r in q['records']],flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--shard',type=int);ap.add_argument('--pilot',action='store_true');a=ap.parse_args()
    assert 1<=a.workers<=4
    if a.pilot:worker(0,1,True);return
    freeze()
    if a.shard is not None:worker(a.shard,a.workers);return
    jobs=[]
    for shard in range(a.workers):
        log=open(ROOT/f'work/integration_{shard}.log','w',encoding='utf-8')
        proc=subprocess.Popen([sys.executable,__file__,'--shard',str(shard),'--workers',str(a.workers)],stdout=log,stderr=subprocess.STDOUT)
        jobs.append((proc,log))
    for proc,log in jobs:
        code=proc.wait();log.close();assert code==0,'Integration worker failed'
    print('Integration bank complete.',flush=True)


if __name__=='__main__':main()
