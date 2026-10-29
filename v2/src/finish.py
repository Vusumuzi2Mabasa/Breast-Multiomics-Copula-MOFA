"""Wait for computation, then perform independent checks and manuscript assembly."""
import argparse,subprocess,sys,time,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true');ap.add_argument('--rscript',required=True);a=ap.parse_args()
    if a.watch:
        while 'computation completed' not in (ROOT/'work/pipeline.log').read_text(errors='replace'):time.sleep(10)
    commands=[('integration_validation',[sys.executable,'v2/src/validate_integration.py']),
              ('integration_R',[a.rscript,'--vanilla','v2/src/check_integration.R','.']),
              ('analysis',[sys.executable,'v2/src/analysis.py']),
              ('figures',[sys.executable,'v2/src/figures.py']),
              ('manuscript',[sys.executable,'v2/src/manuscript.py'])]
    for name,cmd in commands:
        with open(ROOT/f'work/{name}.log','w',encoding='utf-8') as log:
            result=subprocess.run(cmd,cwd=BASE,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise SystemExit(f'{name} failed; inspect v2/work/{name}.log')
        print(name,'passed',flush=True)
    (ROOT/'results/numerical_completion.json').write_text(json.dumps(dict(status='passed',stages=[n for n,c in commands]),indent=2))
