"""Biological integration benchmark with independent seen/unseen-country tests."""
import sys,time
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression,Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score,brier_score_loss,balanced_accuracy_score
from biology import generate
from factors import tune_mofa,transform,fit_preprocess,apply_preprocess,rank_transform_fit,rank_transform_apply,fit_gibbs,project_mofa,copula_projection
sys.path.append(str(Path(__file__).resolve().parents[2]/'src'))
from metrics import full_subspace_recovery,basis

METHODS=['MOFA+ tuned','MOFA+ rank tuned','Gaussian Gibbs matched','Copula Gibbs matched','Rank PCA','Covariates only','Oracle factors']


def subset(d,take):
    n=len(d['y']);out={}
    for key,value in d.items():
        if key=='views':out[key]=[v[take] for v in value]
        elif isinstance(value,np.ndarray) and value.shape and value.shape[0]==n:out[key]=value[take]
        else:out[key]=value
    return out


def datasets(n,seed,scenario,config):
    a=generate(n*2,seed,scenario);train=subset(a,np.flatnonzero(a['country']!=4)[:n])
    b=generate(config['test_size']*2,seed+50000000,scenario);test=subset(b,np.flatnonzero(b['country']!=4)[:config['test_size']])
    c=generate(4000,seed+60000000,scenario);unseen=subset(c,np.flatnonzero(c['country']==4)[:config['unseen_test_size']])
    assert len(train['y'])==n and len(test['y'])==config['test_size'] and len(unseen['y'])==config['unseen_test_size']
    assert 4 not in train['country'] and np.all(unseen['country']==4)
    return train,test,unseen


def classifier(a,b,c,train,test,unseen):
    x=np.c_[train['covariates'],a];xt=np.c_[test['covariates'],b];xu=np.c_[unseen['covariates'],c]
    scaler=StandardScaler().fit(x);fit=LogisticRegression(C=1,max_iter=1000).fit(scaler.transform(x),train['y'])
    p=fit.predict_proba(scaler.transform(xt))[:,1];pu=fit.predict_proba(scaler.transform(xu))[:,1]
    metrics=dict(auc=float(roc_auc_score(test['y'],p)),brier=float(brier_score_loss(test['y'],p)),balanced_accuracy=float(balanced_accuracy_score(test['y'],p>=.5)),unseen_auc=float(roc_auc_score(unseen['y'],pu)),unseen_brier=float(brier_score_loss(unseen['y'],pu)))
    return metrics,p,pu


def protein_error(train_scores,test_scores,train_views,test_views,train,test):
    a=np.c_[train['covariates'],train_scores];b=np.c_[test['covariates'],test_scores]
    scaler=StandardScaler().fit(a);a=scaler.transform(a);b=scaler.transform(b)
    error=[]
    for j in range(16):
        good=np.isfinite(train_views[3][:,j]);testgood=np.isfinite(test_views[3][:,j])
        fit=Ridge(alpha=1).fit(a[good],train_views[3][good,j]);pred=fit.predict(b[testgood])
        error.extend((pred-test_views[3][testgood,j])**2)
    return float(np.mean(error))


def mask_protein(views):
    out=[v.copy() for v in views];out[3][:]=np.nan;return out


def run(train,test,unseen,seed,config):
    records=[];predictions={};models={};tunings={};fixtures={};common_k=None
    raw,params=fit_preprocess(train['views'],train['covariates'])
    rawtest=apply_preprocess(test['views'],test['covariates'],params)
    rawunseen=apply_preprocess(unseen['views'],unseen['covariates'],params)
    rank,margins=rank_transform_fit(raw);ranktest=rank_transform_apply(rawtest,margins);rankunseen=rank_transform_apply(rawunseen,margins)
    matched=[v.copy() for v in raw];matched_test=[v.copy() for v in rawtest];matched_unseen=[v.copy() for v in rawunseen]
    mu=raw[1].mean(0);sd=np.maximum(raw[1].std(0),.15)
    for target,source in [(matched,raw),(matched_test,rawtest),(matched_unseen,rawunseen)]:target[1]=(source[1]-mu)/sd
    for method in METHODS:
        print('Fitting',method,flush=True)
        started=time.perf_counter();rec=dict(method=method,status='ok',diagnostic_pass=True,k=None,rhat_max=None,extended=False,protein_mse=None,capture=None)
        try:
            if method.startswith('MOFA+'):
                isrank='rank' in method;tuned=tune_mofa(train,seed+13,isrank,config)
                model=tuned['model'];k=tuned['k'];va=tuned['train'];vb=transform(tuned,test);vc=transform(tuned,unseen)
                a=project_mofa(model,va,isrank);b=project_mofa(model,vb,isrank);c=project_mofa(model,vc,isrank)
                am=project_mofa(model,mask_protein(va),isrank);bm=project_mofa(model,mask_protein(vb),isrank)
                # Compare protein error on a common observed standardised scale,
                # even when the latent representation used normal scores.
                rec['protein_mse']=protein_error(am,bm,raw,rawtest,train,test)
                rec['k']=k;rec['diagnostic_pass']=model['converged'];tunings[method]=tuned['candidates'];models[method]=tuned
                if not isrank:common_k=k
                fixtures[method.replace(' ','_').replace('+','')+'_scores']=a
            elif 'Gibbs' in method:
                isrank=method.startswith('Copula');assert common_k is not None
                # Gaussian and copula arms use the same raw data, K and priors.
                x=np.concatenate(matched,axis=1)
                model=fit_gibbs(x,common_k,seed+29,isrank,config)
                v,vt,vu=([rank,ranktest,rankunseen] if isrank else [matched,matched_test,matched_unseen])
                a=copula_projection(model['correlation'],np.concatenate(v,axis=1),common_k)
                b=copula_projection(model['correlation'],np.concatenate(vt,axis=1),common_k)
                c=copula_projection(model['correlation'],np.concatenate(vu,axis=1),common_k)
                am=copula_projection(model['correlation'],np.concatenate(mask_protein(v),axis=1),common_k)
                bm=copula_projection(model['correlation'],np.concatenate(mask_protein(vt),axis=1),common_k)
                rec.update(k=common_k,rhat_max=model['rhat_max'],extended=model['extended'],diagnostic_pass=model['converged'])
                rec['protein_mse']=protein_error(am,bm,raw,rawtest,train,test)
                fixtures[method.split()[0]+'_correlation']=model['correlation'];fixtures[method.split()[0]+'_trace']=model['trace']
            elif method=='Rank PCA':
                z=np.nan_to_num(np.concatenate(rank,axis=1));_,_,v=np.linalg.svd(z,full_matrices=False);v=v[:common_k]
                a=z@v.T;b=np.nan_to_num(np.concatenate(ranktest,axis=1))@v.T;c=np.nan_to_num(np.concatenate(rankunseen,axis=1))@v.T
                am=np.nan_to_num(np.concatenate(mask_protein(rank),axis=1))@v.T;bm=np.nan_to_num(np.concatenate(mask_protein(ranktest),axis=1))@v.T
                rec['k']=common_k;rec['protein_mse']=protein_error(am,bm,raw,rawtest,train,test)
            elif method=='Oracle factors':
                a=train['h'];b=test['h'];c=unseen['h'];rec['k']=4
            else:a=np.empty((len(train['y']),0));b=np.empty((len(test['y']),0));c=np.empty((len(unseen['y']),0))
            metrics,p,pu=classifier(a,b,c,train,test,unseen);rec.update(metrics)
            predictions[method]=dict(test=p.tolist(),unseen=pu.tolist())
            if a.shape[1]:rec['capture']=full_subspace_recovery(a,train['h']);rec['rank']=basis(a).shape[1]
        except Exception as exc:
            import traceback
            rec.update(status='failed',error=repr(exc),traceback=traceback.format_exc())
        rec['seconds']=time.perf_counter()-started;records.append(rec)
        print('Finished',method,rec['status'],round(rec['seconds'],2),'seconds',flush=True)
    return dict(records=records,predictions=predictions,tunings=tunings,fixtures=fixtures)
