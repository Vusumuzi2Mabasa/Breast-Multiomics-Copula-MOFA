"""Frozen, resumable simulation runner with independent test participants."""
import argparse, hashlib, json, os, time, traceback
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.linalg import subspace_angles
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from models import fit_copula,fit_mofa,normal_scores_fit,normal_scores_apply
from simulate import generate,standardise,SCENARIOS

ROOT=Path(__file__).resolve().parents[1]


def project_mofa(model,views,gaussian=False):
    # Conditional MAP scores for new participants with fixed fitted loadings.
    n=len(views[0]);k=model['weights'][0].shape[1]
    z=np.zeros((n,k))
    for i in range(n):
        a=np.eye(k);b=np.zeros(k)
        for m in ([0,1,2,3] if gaussian else [0,2,3]):
            good=np.isfinite(views[m][i]);w=model['weights'][m][good]
            tau=np.broadcast_to(np.asarray(model['precision'][m]).reshape(-1),(views[m].shape[1],))[good]
            a+=(w.T*tau)@w;b+=w.T@(tau*views[m][i,good])
        w=model['weights'][1];y=views[1][i]
        zi=np.linalg.solve(a,b)
        if gaussian:z[i]=zi;continue
        for _ in range(15):
            prob=expit(w@zi)
            grad=a@zi-b+w.T@(prob-y)
            hes=a+(w.T*(prob*(1-prob)))@w
            delta=np.linalg.solve(hes,grad);zi-=delta
            if np.max(np.abs(delta))<1e-7:break
        z[i]=zi
    return z


def copula_projection(c,x,k=3):
    vals,vec=np.linalg.eigh(c)
    idx=np.argsort(vals)[-k:]
    load=vec[:,idx]*np.sqrt(np.maximum(vals[idx]-1,.05))
    residual=np.maximum(1-np.sum(load*load,axis=1),.1)
    # Plug-in Gaussian scores with midpoint empirical CDFs for tied categories;
    # missing entries omitted. Exact rank likelihood is used during fitting.
    scores=[]
    for row in x:
        good=np.isfinite(row)
        w=load[good];tau=1/residual[good]
        scores.append(np.linalg.solve(np.eye(k)+(w.T*tau)@w,w.T@(row[good]*tau)))
    return np.asarray(scores)


def predict_metrics(a,b,y,yt,cov,covt):
    x=np.column_stack([cov,a]);xt=np.column_stack([covt,b])
    scale=StandardScaler().fit(x)
    fit=LogisticRegression(C=1,solver='lbfgs',max_iter=500).fit(scale.transform(x),y)
    prob=fit.predict_proba(scale.transform(xt))[:,1]
    return dict(auc=float(roc_auc_score(yt,prob)),ap=float(average_precision_score(yt,prob)),
                brier=float(brier_score_loss(yt,prob)))


def recovery(scores,truth):
    a=scores-scores.mean(0);b=truth-truth.mean(0)
    return float(np.mean(np.cos(subspace_angles(a,b))**2))


def run_one(scenario,n,rep,seed,burn,draws,mofa_iterations,out):
    train=generate(n,seed,scenario);test=generate(200,seed+50000000,scenario)
    va,vb=standardise(train['views'],test['views'])
    x=np.concatenate(train['views'],axis=1);xt=np.concatenate(test['views'],axis=1)
    ns,margins=normal_scores_fit(x);nst=normal_scores_apply(xt,margins)
    # Missing coordinates remain missing for conditional projections.
    ns[~np.isfinite(x)]=np.nan;nst[~np.isfinite(xt)]=np.nan
    records=[];fitoutputs={}
    methods=['MOFA+','MOFA+ Gaussian','Copula','Rank PCA','Clinical only','Oracle factors']
    for method in methods:
        started=time.perf_counter()
        rec=dict(scenario=scenario,n=n,replicate=rep,seed=seed,method=method,
                 status='ok',rhat_max=None,recovery=None,elbo_steps=None)
        try:
            if method.startswith('MOFA+'):
                gaussian=method=='MOFA+ Gaussian'
                av=[v.copy() for v in va];bv=[v.copy() for v in vb]
                if gaussian:
                    mu=av[1].mean(0);sd=np.maximum(av[1].std(0),.1)
                    av[1]=(av[1]-mu)/sd;bv[1]=(bv[1]-mu)/sd
                model=fit_mofa(av,3,seed+13,iterations=mofa_iterations,gaussian=gaussian)
                # Same conditional MAP convention for training and test rows.
                a=project_mofa(model,av,gaussian);b=project_mofa(model,bv,gaussian)
                rec['elbo_steps']=len(model['elbo'])
                rec['diagnostic_pass']=model['converged']
                key='mofa_gaussian' if gaussian else 'mofa'
                fitoutputs[key+'_scores']=a;fitoutputs[key+'_elbo']=model['elbo']
                (out/f'{scenario}_{n}_{rep}_{key}.log').write_text(model['log'],encoding='utf-8')
            elif method=='Copula':
                model=fit_copula(x,3,seed+29,burn=burn,draws=draws)
                rec['extended_sampling']=False
                if not model['converged']:
                    model=fit_copula(x,3,seed+29,burn=burn*3,draws=draws*3)
                    rec['extended_sampling']=True
                a=copula_projection(model['correlation'],ns)
                b=copula_projection(model['correlation'],nst)
                rec['rhat_max']=model['rhat_max'];rec['rhat_median']=model['rhat_median']
                rec['diagnostic_pass']=model['converged']
                fitoutputs['copula_correlation']=model['correlation']
                fitoutputs['copula_trace']=model['trace']
                fitoutputs['copula_scores']=a
            elif method=='Rank PCA':
                fill=np.nan_to_num(ns);_,_,v=np.linalg.svd(fill,full_matrices=False)
                a=fill@v[:3].T;b=np.nan_to_num(nst)@v[:3].T
            elif method=='Oracle factors':a=train['h'];b=test['h']
            else:a=np.empty((n,0));b=np.empty((len(xt),0))
            rec.update(predict_metrics(a,b,train['y'],test['y'],train['covariates'],test['covariates']))
            fitoutputs[method.replace(' ','_').replace('+','')+'_test_metrics']=np.array([rec['auc'],rec['ap'],rec['brier']])
            if a.shape[1]:rec['recovery']=recovery(a,train['h'])
        except Exception as exc:
            rec.update(status='failed',error=repr(exc))
            (out/f'{scenario}_{n}_{rep}_{method.replace("+","")}_error.log').write_text(traceback.format_exc())
        rec['seconds']=time.perf_counter()-started;records.append(rec)
    # Exact generated covariance target applies before technical contamination;
    # no whole-data covariance error claimed for altered/MNAR scenarios.
    fitoutputs.update(y=train['y'],h=train['h'],country=train['country'],data=x,
                      test_y=test['y'],truth=train['truth'])
    np.savez_compressed(out/f'{scenario}_{n}_{rep}_fixture.npz',**fitoutputs)
    return records


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--replicates',type=int)
    ap.add_argument('--shard',type=int,default=0);ap.add_argument('--shards',type=int,default=1)
    args=ap.parse_args();protocol=json.loads((ROOT/'protocol/design.json').read_text())
    folder=ROOT/('work/pilot_verified' if args.pilot else 'results/replicates');folder.mkdir(parents=True,exist_ok=True)
    if not args.pilot:
        lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
        for name,sha in lock['sources'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha, f'Frozen source changed: {name}'
    scenarios=['regular'] if args.pilot else protocol['scenarios']
    ns=[180] if args.pilot else protocol['sample_sizes']
    reps=1 if args.pilot else args.replicates or protocol['replicates']
    for s,scenario in enumerate(scenarios):
        for n in ns:
            for rep in range(reps):
                if rep % args.shards != args.shard:continue
                target=folder/f'{scenario}_{n}_{rep}.json'
                if target.exists():continue
                seed=(protocol['pilot_seed'] if args.pilot else protocol['seed'])+s*100000+n*100+rep
                rows=run_one(scenario,n,rep,seed,protocol['burn'],protocol['draws'],protocol['mofa_iterations'],folder)
                target.write_text(json.dumps(rows,indent=2,allow_nan=False),encoding='utf-8')
                print(scenario,n,rep,[(r['method'],r['status'],round(r.get('auc',0),3),round(r['seconds'],1)) for r in rows],flush=True)


if __name__=='__main__':main()
