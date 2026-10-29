"""Export a small, readable synthetic example independent of evaluation."""
from pathlib import Path
import json
import pandas as pd
import numpy as np
from simulate import generate,NAMES
ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    out=ROOT/'data/example';out.mkdir(parents=True,exist_ok=True)
    d=generate(180,69102902,'missing')
    ids=[f'SIM{i:04d}' for i in range(1,181)]
    for name,view in zip(NAMES,d['views']):
        pd.DataFrame(view,index=ids,columns=[f'{name}_{j+1:03d}' for j in range(view.shape[1])]).rename_axis('synthetic_id').to_csv(out/f'{name}.csv')
    pd.DataFrame(d['raw_counts'],index=ids,columns=[f'RNA_{j+1:03d}' for j in range(24)]).rename_axis('synthetic_id').to_csv(out/'RNA_raw_counts.csv')
    pd.DataFrame(dict(synthetic_id=ids,TNBC=d['y'],country=['C'+str(x+1) for x in d['country']],age_standardised=d['covariates'][:,0],cell_composition_standardised=d['covariates'][:,5])).to_csv(out/'metadata.csv',index=False)
    (out/'manifest.json').write_text(json.dumps(dict(disclosure='Entirely synthetic illustrative data; independent of evaluation and calibration.',seed=d['seed'],participants=180,views=NAMES,missing_protein_fraction=float(np.isnan(d['views'][3]).mean()),genomic_features='Synthetic correlated indicators and dosages; no real variant identifiers or empirical LD.'),indent=2))
    print('Independent example exported to data/example.')
