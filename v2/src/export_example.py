"""Export an independent, small synthetic panel with explicit assay metadata."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from biology import generate
ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    out=ROOT/'data/example';out.mkdir(parents=True,exist_ok=True)
    d=generate(240,901030000,'missing_coordination');ids=[f'SIM_{j+1:04d}' for j in range(240)]
    for key in ['germline','somatic','RNA','protein','counts']:
        a=d[key];pd.DataFrame(a,index=ids,columns=[f'{key}_{j+1:02d}' for j in range(a.shape[1])]).rename_axis('synthetic_id').to_csv(out/f'{key}.csv')
    pd.DataFrame(dict(synthetic_id=ids,TNBC=d['y'],country=['synthetic_C'+str(c+1) for c in d['country']],batch=d['batch'],age_z=d['covariates'][:,0],cell_score=d['covariates'][:,1],treatment=d['covariates'][:,2],whole_assay_library_size=d['library'])).to_csv(out/'metadata.csv',index=False)
    pd.DataFrame([dict(gene=f'SIM_GENE_{j+1:02d}',module=f'MODULE_{j//4+1}',RNA=f'RNA_{j+1:02d}',protein=f'protein_{j+1:02d}',cis_variant=f'germline_{j%8+1:02d}',somatic_feature=f'somatic_{j%8+1:02d}',generating_coordination_beta=float(d['beta'][j])) for j in range(16)]).to_csv(out/'feature_map.csv',index=False)
    (out/'manifest.json').write_text(json.dumps(dict(disclosure='Entirely synthetic, independent of all evaluation and calibration seeds. Feature identifiers do not denote real genes or genomic variants.',participants=240,seed=d['seed'],scenario=d['scenario'],protein_scale='Natural log abundance',RNA_scale='log2 CPM plus one',genomic_LD='Four illustrative latent Gaussian blocks, not empirical African LD',missing_protein_fraction=float(d['missing'].mean())),indent=2))
    print('Independent 240-participant example and gene-module map exported.')
