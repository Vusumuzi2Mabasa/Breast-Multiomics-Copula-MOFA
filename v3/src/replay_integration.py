"""Replay a prespecified complete two-view integration replicate."""
from pathlib import Path
import json
import numpy as np
from simulation import generate
from integration import run
ROOT=Path(__file__).resolve().parents[1]


def main():
    reference=json.loads((ROOT/'results/integration/subtype_eqtl_240_00.json').read_text())
    config=json.loads((ROOT/'protocol/design.json').read_text());seed=reference['seed']
    train=generate(240,seed,'subtype_eqtl');test=generate(config['test_size'],reference['test_seed'],'subtype_eqtl')
    result,_=run(train,test,seed,config)
    differences=[float(np.max(np.abs(np.asarray(p)-np.asarray(reference['predictions'][m])))) for m,p in result['predictions'].items()]
    assert result['selected_k']==reference['selected_k']
    for a,b in zip(result['candidates'],reference['candidates']):assert abs(a['validation_mse']-b['validation_mse'])<1e-10
    assert max(differences)<1e-10
    (ROOT/'results/integration_replay.json').write_text(json.dumps(dict(status='passed',replicate='subtype_eqtl_240_00',max_probability_error=max(differences),methods_checked=len(differences),selected_k=result['selected_k']),indent=2));print('Full integration replay passed')


if __name__=='__main__':main()
