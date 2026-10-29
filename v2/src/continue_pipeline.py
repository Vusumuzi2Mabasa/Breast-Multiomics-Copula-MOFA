"""Bound CPU use: finish dependence first, then integration and official DIABLO."""
import argparse,subprocess,sys,time,datetime,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--wait-dependence',action='store_true');ap.add_argument('--rscript',required=True);ap.add_argument('--workers',type=int,default=4);a=ap.parse_args()
    sources=['v2/src/diablo.R']
    lock=ROOT/'protocol/diablo_lock.json';payload={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in sources}
    if lock.exists():assert json.loads(lock.read_text())['sources']==payload
    else:lock.write_text(json.dumps(dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=payload),indent=2))
    if a.wait_dependence:
        while 'evaluation completed' not in (ROOT/'work/dependence_run.log').read_text(errors='replace'):time.sleep(10)
    else:
        subprocess.run([sys.executable,str(ROOT/'src/run_dependence.py'),'--workers',str(a.workers)],cwd=BASE,check=True)
    with open(ROOT/'work/integration_run.log','w',encoding='utf-8') as log:
        subprocess.run([sys.executable,str(ROOT/'src/run_integration.py'),'--workers',str(a.workers)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    with open(ROOT/'work/diablo_run.log','w',encoding='utf-8') as log:
        subprocess.run([a.rscript,'--vanilla',str(ROOT/'src/diablo.R'),str(BASE)],cwd=BASE,stdout=log,stderr=subprocess.STDOUT,check=True)
    print('Dependence, integration and DIABLO computation completed.',flush=True)
