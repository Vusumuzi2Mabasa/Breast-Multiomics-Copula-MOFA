"""Reproduce Study 3 in a new directory without deleting released results."""
import argparse,subprocess,sys,shutil,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--rscript',default=shutil.which('Rscript'));ap.add_argument('--quarto',default=shutil.which('quarto'));a=ap.parse_args()
    if not a.rscript or not a.quarto:ap.error('Rscript with mixOmics/jsonlite and Quarto are required; supply executable paths if needed.')
    target=a.output.resolve()
    if target.exists():ap.error('Choose a NEW output directory; existing results are never removed.')
    target.mkdir(parents=True)
    for name in ['src','tests','v2/src','v2/tests','v2/protocol','v2/provisional_pca','baseline']:
        if (BASE/name).exists():shutil.copytree(BASE/name,target/name,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['requirements.txt','requirements-lock.txt','LICENSE','CITATION.cff']:shutil.copy2(BASE/name,target/name)
    for name in ['work','v2/work','v2/results','v2/figures','v2/report']:(target/name).mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[name]='1'
    def run(cmd):subprocess.run(cmd,cwd=target,env=env,check=True)
    run([sys.executable,'-m','unittest','discover','-s','tests','-v'])
    run([sys.executable,'-m','unittest','discover','-s','v2/tests','-v'])
    run([sys.executable,'v2/src/run_dependence.py','--workers',str(a.workers)])
    run([sys.executable,'v2/src/run_integration.py','--workers',str(a.workers)])
    run([a.rscript,'--vanilla','v2/src/diablo.R','.'])
    run([sys.executable,'v2/src/validate.py'])
    run([a.rscript,'--vanilla','v2/src/check_conditional.R','.'])
    run([sys.executable,'v2/src/validate_integration.py'])
    run([a.rscript,'--vanilla','v2/src/check_integration.R','.'])
    for name in ['export_example','analysis','figures','manuscript']:run([sys.executable,f'v2/src/{name}.py'])
    run([sys.executable,'v2/src/build_word.py','--quarto',a.quarto])
    run([sys.executable,'v2/src/audit.py','--fresh'])
    print('Fresh scientific reproduction complete:',target)
