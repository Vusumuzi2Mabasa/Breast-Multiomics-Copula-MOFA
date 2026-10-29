"""Independent arithmetic checks and deterministic replay of fixed fixtures."""
import json,hashlib,sys
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss
from simulate import generate,standardise
from models import fit_copula,fit_mofa,normal_scores_fit,normal_scores_apply
from benchmark import project_mofa,copula_projection,recovery
from metrics import full_subspace_recovery
ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/validation';out.mkdir(exist_ok=True)
    p=json.loads((ROOT/'protocol/design.json').read_text())
    checks=[];rfixtures=[]
    for scenario,n in [('regular',180),('missing',360),('heavy_tail',180),('dependence_only',360)]:
        stem=f'{scenario}_{n}_0'
        rows=json.loads((ROOT/f'results/replicates/{stem}.json').read_text())
        seed=rows[0]['seed'];tr=generate(n,seed,scenario);te=generate(200,seed+50000000,scenario)
        f=np.load(ROOT/f'results/replicates/{stem}_fixture.npz')
        np.testing.assert_allclose(np.concatenate(tr['views'],axis=1),f['data'],equal_nan=True)
        np.testing.assert_array_equal(tr['y'],f['y'])
        x=np.concatenate(tr['views'],axis=1);xt=np.concatenate(te['views'],axis=1)
        ns,margin=normal_scores_fit(x);nt=normal_scores_apply(xt,margin)
        ns[~np.isfinite(x)]=np.nan;nt[~np.isfinite(xt)]=np.nan
        va,vb=standardise(tr['views'],te['views'])
        for method in ['Copula','MOFA+']:
            stored=next(r for r in rows if r['method']==method)
            if method=='Copula':
                multiplier=3 if stored['extended_sampling'] else 1
                model=fit_copula(x,3,seed+29,burn=p['burn']*multiplier,draws=p['draws']*multiplier)
                error=float(np.max(np.abs(model['correlation']-f['copula_correlation'])))
                assert error<1e-10,error
                a=copula_projection(model['correlation'],ns);b=copula_projection(model['correlation'],nt)
            else:
                model=fit_mofa(va,3,seed+13,iterations=p['mofa_iterations'])
                a=project_mofa(model,va);b=project_mofa(model,vb)
                error=float(np.max(np.abs(a-f['mofa_scores'])))
                assert error<1e-8,error
            scale=StandardScaler().fit(np.c_[tr['covariates'],a])
            logit=LogisticRegression(C=1,solver='lbfgs',max_iter=500).fit(scale.transform(np.c_[tr['covariates'],a]),tr['y'])
            prob=logit.predict_proba(scale.transform(np.c_[te['covariates'],b]))[:,1]
            metrics=[roc_auc_score(te['y'],prob),average_precision_score(te['y'],prob),brier_score_loss(te['y'],prob),recovery(a,tr['h'])]
            maxerror=max(abs(metrics[j]-stored[m]) for j,m in enumerate(['auc','ap','brier','recovery']))
            assert maxerror<1e-9,maxerror
            checks.append(dict(fixture=stem,method=method,model_error=error,metric_error=maxerror))
            metrics[-1]=full_subspace_recovery(a,tr['h'])
            rfixtures.append(dict(name=stem+'_'+method,truth=te['y'].tolist(),probability=prob.tolist(),scores=a.tolist(),latent=tr['h'].tolist(),expected=metrics))
    (out/'replay.json').write_text(json.dumps(dict(status='passed',checks=checks),indent=2))
    (out/'R_fixtures.json').write_text(json.dumps(rfixtures))
    print('Eight deterministic model and metric replays passed.')


if __name__=='__main__':main()
