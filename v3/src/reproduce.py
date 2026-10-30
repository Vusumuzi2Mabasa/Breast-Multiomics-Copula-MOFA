"""Reproduce the genotype/RNA experiment in a NEW directory."""
from pathlib import Path
import argparse,subprocess,shutil,sys,os
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--rscript',default=shutil.which('Rscript'));ap.add_argument('--quarto',default=shutil.which('quarto'));a=ap.parse_args()
    if not a.rscript or not a.quarto:ap.error('Rscript with mixOmics/jsonlite and Quarto required')
    target=a.output.resolve()
    if target.exists():ap.error('Choose a new output directory; released results are never deleted')
    target.mkdir(parents=True)
    for name in ['v3/src','v3/tests','v3/protocol','v3/baseline']:shutil.copytree(BASE/name,target/name,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['requirements.txt','requirements-lock.txt','LICENSE']:shutil.copy2(BASE/name,target/name)
    for name in ['work','v3/work','v3/results','v3/figures','v3/report','v3/data/example']:(target/name).mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[key]='1'
    def run(cmd):subprocess.run(cmd,cwd=target,env=env,check=True)
    run([sys.executable,'-m','unittest','discover','-s','v3/tests','-v'])
    for stage in ['eqtl','integration']:run([sys.executable,'v3/src/run.py','--stage',stage,'--workers',str(a.workers)])
    run([a.rscript,'--vanilla','v3/src/diablo.R','.'])
    run([sys.executable,'v3/src/replay_integration.py'])
    run([sys.executable,'v3/src/analysis.py']);run([a.rscript,'--vanilla','v3/src/check_R.R','.'])
    run([sys.executable,'v3/src/robustness.py']);run([a.rscript,'--vanilla','v3/src/check_HC3.R','.'])
    for name in ['figures','manuscript']:run([sys.executable,f'v3/src/{name}.py'])
    run([sys.executable,'v3/src/build_word.py','--quarto',a.quarto]);run([sys.executable,'v3/src/audit.py','--fresh'])
    print('Fresh scientific reproduction complete:',target)


if __name__=='__main__':main()
