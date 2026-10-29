"""Create a fresh output directory and rerun the complete benchmark safely."""
import argparse,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--workers',type=int,default=2)
    ap.add_argument('--rscript',default=shutil.which('Rscript'))
    ap.add_argument('--quarto',default=shutil.which('quarto'))
    ap.add_argument('--numerical-only',action='store_true')
    a=ap.parse_args()
    if not a.numerical_only and (not a.rscript or not a.quarto):
        ap.error('Full reproduction requires Rscript and Quarto; set paths or request --numerical-only.')
    dest=a.output.resolve()
    if dest.exists():ap.error('Use a new output directory; existing files are never removed.')
    dest.mkdir(parents=True)
    for name in ['src','tests','protocol']:
        shutil.copytree(ROOT/name,dest/name,ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['requirements.txt','environment.yml','LICENSE','CITATION.cff']:
        shutil.copy2(ROOT/name,dest/name)
    for name in ['results','figures','report','work']:(dest/name).mkdir()
    env=os.environ.copy()
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:env[name]='1'
    def run(args):subprocess.run(args,cwd=dest,env=env,check=True)
    run([sys.executable,'-m','unittest','discover','-s','tests','-v'])
    run([sys.executable,'src/run_all.py','--workers',str(a.workers)])
    run([sys.executable,'src/dependence_validation.py'])
    run([sys.executable,'src/export_example.py'])
    run([sys.executable,'src/analyse.py'])
    run([sys.executable,'src/validate.py'])
    if a.rscript:run([a.rscript,'--vanilla','src/check_metrics.R','.'])
    run([sys.executable,'src/manuscript.py'])
    if not a.numerical_only:run([sys.executable,'src/build_word.py','--quarto',a.quarto])
    print('Fresh reproduction complete:',dest)
