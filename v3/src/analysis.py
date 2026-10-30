"""Aggregate every prespecified replicate, retaining diagnostics and failures."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score,roc_auc_score,brier_score_loss
from simulation import generate
from eqtl import analyse
ROOT=Path(__file__).resolve().parents[1]


def meanse(values):
    a=np.asarray(values,float);a=a[np.isfinite(a)]
    return (float(a.mean()),float(a.std(ddof=1)/np.sqrt(len(a)))) if len(a)>1 else (float(a[0]) if len(a) else None,None)


def main():
    eq=[];gene=[];failures=[]
    for p in sorted((ROOT/'results/eqtl').glob('*.json')):
        d=json.loads(p.read_text())
        if d['status']!='ok':failures.append(d);continue
        for family,truth in [('omnibus','true_eqtl'),('interaction','true_interaction')]:
            selected=[r[family+'_q']<.05 for r in d['genes']];t=np.array([r[truth] for r in d['genes']]);s=np.array(selected)
            eq.append(dict(scenario=d['scenario'],n=d['n'],replicate=d['replicate'],seed=d['seed'],family=family,discoveries=int(s.sum()),true_discoveries=int((s&t).sum()),false_discoveries=int((s&~t).sum()),sensitivity=float((s&t).sum()/t.sum()) if t.sum() else np.nan,FDP=float((s&~t).sum()/max(s.sum(),1)),any_false=float((s&~t).any()),TNBC=d['TNBC']))
        for row in d['genes']:gene.append(dict(scenario=d['scenario'],n=d['n'],replicate=d['replicate'],**row))
    q=pd.DataFrame(eq);q.to_csv(ROOT/'results/eqtl_evaluation.csv',index=False);pd.DataFrame(gene).to_csv(ROOT/'results/eqtl_gene_records.csv',index=False)
    summaries=[]
    for key,g in q.groupby(['scenario','n','family']):
        row=dict(zip(['scenario','n','family'],key));row['replicates']=len(g)
        for metric in ['discoveries','sensitivity','FDP','any_false']:
            mu,se=meanse(g[metric]);row[metric]=mu;row[metric+'_mcse']=se
        summaries.append(row)
    pd.DataFrame(summaries).to_csv(ROOT/'results/eqtl_summary.csv',index=False)
    records=[];fixtures=[]
    for p in sorted((ROOT/'results/integration').glob('*.json')):
        d=json.loads(p.read_text())
        if d['status']!='ok':failures.append(d);continue
        for row in d['records']:
            r=dict(scenario=d['scenario'],n=d['n'],replicate=d['replicate'],seed=d['seed'],**row);records.append(r)
            y=np.array(d['test_y']);prob=np.array(d['predictions'][row['method']]);assert len(y)==len(prob)==360
            assert abs(roc_auc_score(y,prob)-row['auc'])<1e-12
            assert abs(balanced_accuracy_score(y,prob>=.5)-row['balanced_accuracy'])<1e-12
            assert abs(brier_score_loss(y,prob)-row['brier'])<1e-12
            if d['replicate']==0:fixtures.append(dict(method=row['method'],y=y.tolist(),p=prob.tolist(),auc=row['auc'],balanced_accuracy=row['balanced_accuracy'],brier=row['brier']))
        for c in d['candidates']:
            assert len({x['hash'] for x in c['initializations']})==2
        assert d['selected_k']==min(d['candidates'],key=lambda c:(c['validation_mse'],c['k']))['k']
    for p in sorted((ROOT/'results/diablo').glob('*.json')):
        d=json.loads(p.read_text())
        if d['status']!='ok':failures.append(d);continue
        original=json.loads((ROOT/f'results/integration/{p.name}').read_text());assert d['test_y']==original['test_y']
        assert abs(balanced_accuracy_score(d['test_y'],d['test_class'])-d['balanced_accuracy'])<1e-12
        records.append(dict(scenario=d['scenario'],n=d['n'],replicate=d['replicate'],seed=d['seed'],method='DIABLO',status='ok',balanced_accuracy=d['balanced_accuracy'],diagnostic_pass=d['diagnostic_pass'],k=d['k'],auc=np.nan,brier=np.nan))
    df=pd.DataFrame(records);df.to_csv(ROOT/'results/integration_evaluation.csv',index=False);rows=[]
    for key,g in df.groupby(['scenario','n','method']):
        row=dict(zip(['scenario','n','method'],key));row.update(replicates=len(g),diagnostic_flags=int((~g.diagnostic_pass).sum()))
        for metric in ['balanced_accuracy','auc','brier']:
            row[metric],row[metric+'_mcse']=meanse(g[metric])
        rows.append(row)
    pd.DataFrame(rows).to_csv(ROOT/'results/integration_summary.csv',index=False)
    (ROOT/'results/R_metric_fixtures.json').write_text(json.dumps(fixtures))
    # Eight R likelihood fixtures span all scenarios and sample sizes; seeds never selected by results.
    qfixtures=[]
    for n in [240,480]:
        for sc in ['null','shared_eqtl','subtype_eqtl','abundance_only']:
            meta=json.loads((ROOT/f'results/eqtl/{sc}_{n}_000.json').read_text());d=generate(n,meta['seed'],sc)
            qfixtures.append(dict(scenario=sc,n=n,seed=meta['seed'],genotype=d['genotype'].tolist(),RNA=d['RNA'].tolist(),covariates=d['covariates'].tolist(),y=d['y'].tolist(),expected=meta['genes']))
    (ROOT/'results/R_eqtl_fixtures.json').write_text(json.dumps(qfixtures))
    # Replay all genes from eight complete cohorts.
    for d in qfixtures:
        current=analyse(generate(d['n'],d['seed'],d['scenario']))
        assert current==d['expected']
    flags=int((~df.diagnostic_pass).sum())
    audit=dict(status='passed' if not failures else 'failed',eqtl_cohorts=len(q)//2,eqtl_gene_records=len(gene),integration_cohorts=len(list((ROOT/'results/integration').glob('*.json'))),integration_method_records=len(df),diagnostic_flags=flags,failures=failures,replayed_eqtl_cohorts=8,checked_probability_records=len(fixtures))
    (ROOT/'results/numerical_audit.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2));assert not failures


if __name__=='__main__':main()
