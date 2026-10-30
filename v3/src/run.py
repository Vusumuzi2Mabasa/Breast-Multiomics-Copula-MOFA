"""Resumable independent simulations; scientific sources locked before evaluation."""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys,os,time,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from simulation import generate,SCENARIOS,export_example
from eqtl import analyse
from integration import run as integrate
ROOT=Path(__file__).resolve().parents[1]


def freeze():
    paths=['protocol/design.json','src/simulation.py','src/eqtl.py','src/integration.py','src/run.py','src/diablo.R']
    content={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    lock=ROOT/'protocol/source_lock.json'
    if lock.exists():assert json.loads(lock.read_text())['sha256']==content,'Scientific sources changed after source lock'
    else:lock.write_text(json.dumps(dict(sha256=content),indent=2))


def eqtl_job(spec):
    scenario,n,rep=spec;seed=603000000+SCENARIOS.index(scenario)*1000000+n*1000+rep
    out=ROOT/f'results/eqtl/{scenario}_{n}_{rep:03}.json'
    if out.exists():return
    try:
        d=generate(n,seed,scenario);rows=analyse(d)
        result=dict(status='ok',scenario=scenario,n=n,replicate=rep,seed=seed,TNBC=int(d['y'].sum()),genes=rows)
    except Exception as exc:result=dict(status='failed',scenario=scenario,n=n,replicate=rep,seed=seed,error=repr(exc),traceback=traceback.format_exc())
    out.write_text(json.dumps(result,indent=2))


def integration_job(spec):
    scenario,n,rep=spec;seed=803000000+SCENARIOS.index(scenario)*1000000+n*1000+rep
    out=ROOT/f'results/integration/{scenario}_{n}_{rep:02}.json'
    if out.exists():return
    config=json.loads((ROOT/'protocol/design.json').read_text())
    try:
        train=generate(n,seed,scenario);test=generate(config['test_size'],seed+100000000,scenario)
        result,fixture=integrate(train,test,seed,config)
        result.update(status='ok',scenario=scenario,n=n,replicate=rep,seed=seed,test_seed=seed+100000000,test_y=test['y'].tolist(),train_TNBC=int(train['y'].sum()))
        if rep==0:
            np.savez_compressed(out.with_suffix('.npz'),**fixture,train_y=train['y'],test_y=test['y'],h=train['h'])
        folder=ROOT/f'work/diablo_inputs/{out.stem}';folder.mkdir(parents=True,exist_ok=True)
        for label,d in [('train',train),('test',test)]:
            for key in ['genotype','RNA','covariates','y']:np.savetxt(folder/f'{label}_{key}.csv',d[key],delimiter=',',fmt='%.16g')
    except Exception as exc:result=dict(status='failed',scenario=scenario,n=n,replicate=rep,seed=seed,error=repr(exc),traceback=traceback.format_exc())
    out.write_text(json.dumps(result,indent=2))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['pilot','eqtl','integration'],required=True);ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
    for folder in ['results/eqtl','results/integration','results/diablo','work','figures','report']:(ROOT/folder).mkdir(parents=True,exist_ok=True)
    config=json.loads((ROOT/'protocol/design.json').read_text())
    if a.stage=='pilot':
        # Dedicated development seed; never included in the evaluation banks.
        d=generate(120,991030001,'subtype_eqtl');test=generate(120,991030002,'subtype_eqtl')
        rows=analyse(d);result,_=integrate(d,test,991030001,config)
        (ROOT/'work/pilot.json').write_text(json.dumps(dict(qtl=rows,integration=result),indent=2))
        meta=ROOT/'work/pilot_integration';meta.mkdir(exist_ok=True)
        (meta/'pilot.json').write_text(json.dumps(dict(status='ok',seed=991030001,scenario='pilot',n=120,replicate=0)))
        folder=ROOT/'work/pilot_diablo_inputs/pilot';folder.mkdir(parents=True,exist_ok=True)
        for label,data in [('train',d),('test',test)]:
            for key in ['genotype','RNA','covariates','y']:np.savetxt(folder/f'{label}_{key}.csv',data[key],delimiter=',',fmt='%.16g')
        print('Pilot completed',flush=True);return
    freeze();reps=config['eqtl_replicates' if a.stage=='eqtl' else 'integration_replicates']
    specs=[(s,n,r) for s in config['scenarios'] for n in config['sample_sizes'] for r in range(reps)]
    job=eqtl_job if a.stage=='eqtl' else integration_job
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(job,s) for s in specs]
        for k,future in enumerate(as_completed(futures),1):future.result();print(a.stage,k,'/',len(specs),flush=True)
    export_example();print(a.stage,'complete',flush=True)


if __name__=='__main__':main()
