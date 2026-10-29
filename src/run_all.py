"""Run independent simulation shards with bounded CPU use and resumable output."""
import argparse,subprocess,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=2);args=ap.parse_args()
    assert 1<=args.workers<=4
    (ROOT/'work').mkdir(exist_ok=True)
    jobs=[]
    for shard in range(args.workers):
        log=open(ROOT/f'work/shard_{shard}.log','w',encoding='utf-8')
        job=subprocess.Popen([sys.executable,str(ROOT/'src/benchmark.py'),'--shard',str(shard),'--shards',str(args.workers)],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        jobs.append((job,log))
    codes=[]
    for job,log in jobs:codes.append(job.wait());log.close()
    print('Shard exit codes:',codes,flush=True)
    if any(codes):raise SystemExit(1)
    subprocess.run([sys.executable,str(ROOT/'src/dependence.py')],cwd=ROOT,check=True)
