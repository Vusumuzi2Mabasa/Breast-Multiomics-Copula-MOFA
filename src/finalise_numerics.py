"""Wait for the fixed evaluation bank, then validate and assemble its outputs."""
import argparse,json,subprocess,sys,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true')
    ap.add_argument('--rscript',default=shutil.which('Rscript'));a=ap.parse_args()
    if not a.rscript:ap.error('Rscript required for independent verification.')
    if a.watch:
        while len(list((ROOT/'results/replicates').glob('*.json')))<360:
            time.sleep(10)
        while 'Finished 1000 calibration' not in (ROOT/'work/run_all.log').read_text(errors='replace'):
            time.sleep(5)
    steps=[('analysis',[sys.executable,'src/analyse.py']),
           ('validation',[sys.executable,'src/validate.py']),
           ('R_check',[a.rscript,'--vanilla','src/check_metrics.R','.']),
           ('manuscript',[sys.executable,'src/manuscript.py'])]
    for name,cmd in steps:
        with open(ROOT/f'work/{name}.log','w',encoding='utf-8') as log:
            result=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise SystemExit(f'{name} failed; inspect work/{name}.log')
        print(name,'passed',flush=True)
    (ROOT/'results/numerical_completion.json').write_text(json.dumps(dict(status='passed',stages=[s[0] for s in steps]),indent=2))
