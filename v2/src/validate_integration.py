"""Audit tuned models, independent-test arithmetic, DIABLO and one full replay."""
import json,sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score,brier_score_loss,balanced_accuracy_score
from integration import datasets,run
ROOT=Path(__file__).resolve().parents[1]


def main():
    out=ROOT/'results/validation';out.mkdir(exist_ok=True)
    fixtures=[];checks=[]
    for path in sorted((ROOT/'results/integration').glob('*.json')):
        d=json.loads(path.read_text());assert len(d['records'])==7
        k={r['method']:r['k'] for r in d['records']}
        assert k['Gaussian Gibbs matched']==k['Copula Gibbs matched']==k['Rank PCA']==k['MOFA+ tuned']
        for method,candidates in d['tunings'].items():
            assert len(candidates)==4
            for candidate in candidates:
                assert len({s['sha256'] for s in candidate['initializations']})==2
                assert abs(candidate['elbo']-max(s['elbo'] for s in candidate['initializations']))<1e-10
            selected=min(candidates,key=lambda c:(c['validation_mse'],c['k']))['k'];assert k[method]==selected
        y=np.array(d['test_y']);yu=np.array(d['unseen_y'])
        for r in d['records']:
            assert r['status']=='ok',r
            p=np.array(d['predictions'][r['method']]['test']);pu=np.array(d['predictions'][r['method']]['unseen'])
            error=max(abs(roc_auc_score(y,p)-r['auc']),abs(brier_score_loss(y,p)-r['brier']),abs(roc_auc_score(yu,pu)-r['unseen_auc']),abs(balanced_accuracy_score(y,p>=.5)-r['balanced_accuracy']))
            assert error<1e-12
            if d['replicate']==0:fixtures.append(dict(name=path.stem+'_'+r['method'],truth=y.tolist(),p=p.tolist(),unseen_truth=yu.tolist(),unseen_p=pu.tolist(),expected=[r['auc'],r['brier'],r['unseen_auc'],r['balanced_accuracy']]))
        di=json.loads((ROOT/'results/diablo'/path.name).read_text());assert di['status']=='ok',di
        assert abs(balanced_accuracy_score(y,di['test_class'])-di['balanced_accuracy'])<1e-12
        assert abs(balanced_accuracy_score(yu,di['unseen_class'])-di['unseen_balanced_accuracy'])<1e-12
        checks.append(dict(fixture=path.stem,metrics='passed',initializations='distinct',factor_budget='matched',diablo='passed'))
    assert len(checks)==96
    # Replay one prespecified biological integration cohort, including tuning.
    p=json.loads((ROOT/'protocol/design.json').read_text())['integration'];p['factor_grid']=[2,4,6,10]
    stored=json.loads((ROOT/'results/integration/both_240_00.json').read_text());data=datasets(240,stored['seed'],'both',p)
    q=run(*data,stored['seed'],p);errors=[]
    for a,b in zip(q['records'],stored['records']):
        assert a['method']==b['method'] and a['status']=='ok'
        errors.append(max(abs(a[m]-b[m]) for m in ['auc','brier','balanced_accuracy','unseen_auc']))
    assert max(errors)<1e-9,max(errors)
    (out/'integration_checks.json').write_text(json.dumps(dict(status='passed',replicates=96,full_replays=1,max_replay_error=max(errors),checks=checks),indent=2))
    (out/'R_integration_fixtures.json').write_text(json.dumps(fixtures))
    print('Integration selection, paired-test arithmetic, 96 DIABLO results and full replay passed.')


if __name__=='__main__':main()
