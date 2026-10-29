"""Aggregate all evaluation records without excluding difficult settings."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from dependence import METHODS
ROOT=Path(__file__).resolve().parents[1]


def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    centre=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(centre-half),float(centre+half)


def load_dep(phase):
    groups={};failures=[]
    for path in sorted((ROOT/f'results/dependence_chunks/{phase}').glob('*.npz')):
        meta=json.loads(path.with_suffix('.json').read_text())
        assert hashlib.sha256(path.read_bytes()).hexdigest()==meta['sha256']
        data=np.load(path);key=(meta['n'],meta['scenario']);groups.setdefault(key,[]).append({name:data[name] for name in data.files})
        failures.extend([dict(n=meta['n'],scenario=meta['scenario'],**v) for v in meta['errors']])
    return {key:{name:np.concatenate([c[name] for c in chunks],axis=0) for name in chunks[0]} for key,chunks in groups.items()},failures


def outcomes(stat,threshold,active):
    # Failed evaluations remain nondiscoveries; failures are reported separately.
    q=np.nan_to_num(stat,nan=0)**2
    gene=q>np.array(threshold['gene'])[None,:,None]
    module=q.reshape(len(q),len(METHODS),4,4).sum(3)>np.array(threshold['module'])[None,:,None]
    if active:
        return dict(module_detection=module[:,:,0],false_module=module[:,:,1:].any(2),gene_sensitivity=gene[:,:,:4].mean(2),false_gene=gene[:,:,4:].any(2))
    return dict(module_detection=np.full(module.shape[:2],np.nan),false_module=module.any(2),gene_sensitivity=np.full(gene.shape[:2],np.nan),false_gene=gene.any(2))


def dependence_summary():
    design=json.loads((ROOT/'protocol/design.json').read_text());p=design['dependence']
    banks,calfail=load_dep('calibration');evaluation,evalfail=load_dep('evaluation');assert not calfail
    thresholds=json.loads((ROOT/'results/family_thresholds.json').read_text())['thresholds']
    summary=[];records=[]
    seeds=[]
    for (n,scenario),g in evaluation.items():
        seeds.extend(g['seeds'].tolist())
        active=scenario in p['alternative_scenarios'] or scenario=='heavy_tail_coordination'
        o=outcomes(g['statistics'],thresholds[str(n)],active)
        for mi,method in enumerate(METHODS):
            row=dict(n=n,scenario=scenario,method=method,evaluations=len(g['statistics']),failures=int(np.isnan(g['statistics'][:,mi]).any(1).sum()),mean_seconds_per_dataset=float(np.nanmean(g['seconds'])))
            for metric,a in o.items():
                values=a[:,mi]
                row[metric]=float(np.mean(values)) if np.isfinite(values).any() else np.nan
                if metric!='gene_sensitivity' and np.isfinite(values).any():row[metric+'_lower'],row[metric+'_upper']=wilson(sum(values),len(values))
            summary.append(row)
        for i,seed in enumerate(g['seeds']):
            for mi,method in enumerate(METHODS):
                records.append(dict(n=n,scenario=scenario,seed=int(seed),method=method,**{metric:float(a[i,mi]) for metric,a in o.items()}))
    assert len(seeds)==len(set(seeds))
    pd.DataFrame(summary).to_csv(ROOT/'results/dependence_summary.csv',index=False)
    pd.DataFrame(records).to_csv(ROOT/'results/dependence_evaluation.csv',index=False)
    # Paired, calibration-aware contrasts; indices are shared by methods.
    rng=np.random.default_rng(501030000);contrast=[];draws=[]
    for n in design['sample_sizes']:
        bank=[banks[(n,s)]['statistics'] for s in p['null_scenarios']]
        for scenario in p['alternative_scenarios']:
            ev=evaluation[(n,scenario)]['statistics'];point=outcomes(ev,thresholds[str(n)],True)['module_detection']
            boot=np.empty((250,len(METHODS)))
            for b in range(250):
                cut=[]
                for stat in bank:
                    q=stat[rng.integers(0,len(stat),len(stat))]**2
                    family=q.reshape(len(q),len(METHODS),4,4).sum(3).max(2)
                    cut.append(np.sort(family,axis=0)[-25])
                cut=np.max(cut,axis=0)
                q=ev[rng.integers(0,len(ev),len(ev))]**2
                detected=q.reshape(len(q),len(METHODS),4,4).sum(3)[:,:,0]>cut[None,:]
                boot[b]=detected.mean(0)
            for comp in ['Product regression','Normal-score difference','Partial Spearman']:
                a=METHODS.index('Conditional copula');c=METHODS.index(comp)
                delta=boot[:,a]-boot[:,c];lo,hi=np.quantile(delta,[.025,.975])
                contrast.append(dict(n=n,scenario=scenario,comparator=comp,difference=float(np.mean(point[:,a].astype(float)-point[:,c].astype(float))),lower=float(lo),upper=float(hi),bootstrap_replicates=250))
                draws.extend(dict(n=n,scenario=scenario,comparator=comp,draw=i,difference=float(v)) for i,v in enumerate(delta))
    pd.DataFrame(contrast).to_csv(ROOT/'results/dependence_contrasts.csv',index=False)
    pd.DataFrame(draws).to_csv(ROOT/'results/dependence_bootstrap.csv',index=False)
    audit=dict(calibration_datasets=sum(len(v['seeds']) for v in banks.values()),evaluation_datasets=len(seeds),calibration_failures=calfail,evaluation_failures=evalfail,methods=METHODS,
               total_method_records=len(records),threshold_sha256=hashlib.sha256((ROOT/'results/family_thresholds.json').read_bytes()).hexdigest())
    (ROOT/'results/dependence_audit.json').write_text(json.dumps(audit,indent=2))
    print('Dependence aggregation:',audit['calibration_datasets'],audit['evaluation_datasets'],'datasets',flush=True)


def integration_summary():
    rows=[];tuning=[]
    for path in sorted((ROOT/'results/integration').glob('*.json')):
        q=json.loads(path.read_text())
        for r in q['records']:
            if r['method']=='Covariates only' and r['status']=='ok':
                from integration import protein_error
                from factors import fit_preprocess,apply_preprocess
                folder=ROOT/'work/diablo_inputs'/path.stem
                x=np.loadtxt(folder/'train_covariates.csv',delimiter=',');xt=np.loadtxt(folder/'test_covariates.csv',delimiter=',')
                a=[np.loadtxt(folder/f'train_{m}.csv',delimiter=',') for m in ['germline','somatic','RNA','protein']]
                b=[np.loadtxt(folder/f'test_{m}.csv',delimiter=',') for m in ['germline','somatic','RNA','protein']]
                av,params=fit_preprocess(a,x);bv=apply_preprocess(b,xt,params)
                r['protein_mse']=protein_error(np.empty((len(x),0)),np.empty((len(xt),0)),av,bv,dict(covariates=x),dict(covariates=xt))
            rows.append(dict(scenario=q['scenario'],n=q['n'],replicate=q['replicate'],seed=q['seed'],**r))
        for method,items in q['tunings'].items():
            for t in items:tuning.append(dict(scenario=q['scenario'],n=q['n'],replicate=q['replicate'],method=method,**t))
    df=pd.DataFrame(rows);assert len(df)==96*7,len(df)
    df.to_csv(ROOT/'results/integration_evaluation.csv',index=False);pd.DataFrame(tuning).to_csv(ROOT/'results/integration_tuning.csv',index=False)
    metric=['auc','brier','balanced_accuracy','unseen_auc','unseen_brier','capture','protein_mse','seconds']
    summary=[]
    for (s,n,m),g in df.groupby(['scenario','n','method']):
        good=g[g.status=='ok'];r=dict(scenario=s,n=n,method=m,planned=len(g),successful=len(good),diagnostic_flags=int((good.diagnostic_pass==False).sum()))
        for name in metric:
            a=good[name].dropna().to_numpy();r[name]=float(np.mean(a)) if len(a) else np.nan;r[name+'_mcse']=float(np.std(a,ddof=1)/np.sqrt(len(a))) if len(a)>1 else np.nan
        summary.append(r)
    pd.DataFrame(summary).to_csv(ROOT/'results/integration_summary.csv',index=False)
    rng=np.random.default_rng(601030000);contrasts=[];sensitivity=[]
    for (s,n),g in df[df.status=='ok'].groupby(['scenario','n']):
        for comp in ['MOFA+ tuned','MOFA+ rank tuned','Gaussian Gibbs matched','Rank PCA']:
            for name in ['auc','capture','protein_mse','unseen_auc']:
                w=g.pivot(index='replicate',columns='method',values=name).dropna(subset=['Copula Gibbs matched',comp])
                a=(w['Copula Gibbs matched']-w[comp]).to_numpy()*(-1 if name=='protein_mse' else 1)
                b=a[rng.integers(0,len(a),(2000,len(a)))].mean(1);lo,hi=np.quantile(b,[.025,.975])
                contrasts.append(dict(scenario=s,n=n,metric=name,comparator=comp,difference=np.mean(a),lower=lo,upper=hi,n_pairs=len(a)))
                diagnostics=g.pivot(index='replicate',columns='method',values='diagnostic_pass')
                good=diagnostics['Copula Gibbs matched']&diagnostics[comp]
                selected=w.loc[w.index.intersection(good[good].index)]
                delta=(selected['Copula Gibbs matched']-selected[comp]).to_numpy()*(-1 if name=='protein_mse' else 1)
                sensitivity.append(dict(scenario=s,n=n,metric=name,comparator=comp,all_pairs=len(a),diagnostic_pairs=len(delta),all_difference=float(np.mean(a)),diagnostic_difference=float(np.mean(delta)) if len(delta) else np.nan))
    pd.DataFrame(contrasts).to_csv(ROOT/'results/integration_contrasts.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(ROOT/'results/diagnostic_sensitivity.csv',index=False)
    transport=[]
    for (s,n,m),g in df[df.status=='ok'].groupby(['scenario','n','method']):
        delta=(g.unseen_auc-g.auc).to_numpy();draw=delta[rng.integers(0,len(delta),(2000,len(delta)))].mean(1)
        lo,hi=np.quantile(draw,[.025,.975])
        transport.append(dict(scenario=s,n=n,method=m,difference=np.mean(delta),lower=lo,upper=hi,replicates=len(delta)))
    pd.DataFrame(transport).to_csv(ROOT/'results/transport_summary.csv',index=False)
    diablo=[json.loads(p.read_text()) for p in sorted((ROOT/'results/diablo').glob('*.json'))];assert len(diablo)==96
    pd.DataFrame([{k:v for k,v in d.items() if k not in ['candidates','test_class','unseen_class','software','warnings']} for d in diablo]).to_csv(ROOT/'results/diablo_evaluation.csv',index=False)
    (ROOT/'results/integration_audit.json').write_text(json.dumps(dict(training_test_pairs=df.seed.nunique(),method_records=len(df),fit_failures=int((df.status!='ok').sum()),diagnostic_flags=int((df.diagnostic_pass==False).sum()),diablo_fits=len(diablo),diablo_failures=sum(d['status']!='ok' for d in diablo)),indent=2))
    print('Integration aggregation:',df.seed.nunique(),'pairs',flush=True)


if __name__=='__main__':
    if '--dependence-only' in sys.argv:dependence_summary()
    else:dependence_summary();integration_summary()
