"""Fairer integration: technical adjustment, train-only tuning, matched samplers."""
import contextlib,io,sys,time,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from scipy.special import ndtri
sys.path.append(str(Path(__file__).resolve().parents[2]/'src'))
from models import rank_metadata,pack_metadata,latent_update_packed,assert_order,correlation,split_rhat,normal_scores_fit,normal_scores_apply
from benchmark import project_mofa,copula_projection


def fit_preprocess(views,covariates):
    x=np.c_[np.ones(len(covariates)),covariates];params=[];out=[]
    for j,a in enumerate(views):
        if j<2:
            b=np.zeros((x.shape[1],a.shape[1]));b[0]=np.nanmean(a,axis=0)
        else:
            b=np.zeros((x.shape[1],a.shape[1]))
            for c in range(a.shape[1]):
                ok=np.isfinite(a[:,c]);b[:,c]=np.linalg.lstsq(x[ok],a[ok,c],rcond=None)[0]
        residual=a-x@b
        sd=np.maximum(np.nanstd(residual,axis=0),.15)
        out.append((residual/sd) if j!=1 else a.copy())
        params.append((b,sd))
    return out,params


def apply_preprocess(views,covariates,params):
    x=np.c_[np.ones(len(covariates)),covariates]
    return [(a-x@b)/sd if j!=1 else a.copy() for j,(a,(b,sd)) in enumerate(zip(views,params))]


def rank_transform_fit(views):
    flat=np.concatenate(views,axis=1);z,margins=normal_scores_fit(flat)
    z[~np.isfinite(flat)]=np.nan
    return list(np.split(z,np.cumsum([8,8,16]),axis=1)),margins


def rank_transform_apply(views,margins):
    flat=np.concatenate(views,axis=1);z=normal_scores_apply(flat,margins)
    z[~np.isfinite(flat)]=np.nan
    return list(np.split(z,np.cumsum([8,8,16]),axis=1))


def official_mofa(views,k,seed,rank=False,iterations=500,perturb=False):
    from mofapy2.run.entry_point import entry_point
    log=io.StringIO();start=time.perf_counter()
    with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        e=entry_point();e.set_data_options(scale_views=False,center_groups=False)
        e.set_data_matrix([[v.copy()] for v in views],likelihoods=['gaussian']*4 if rank else ['gaussian','bernoulli','gaussian','gaussian'],views_names=['germline','somatic','RNA','protein'])
        e.set_model_options(factors=k,ard_weights=True,spikeslab_weights=True,ard_factors=False,spikeslab_factors=False)
        e.set_train_options(iter=iterations,convergence_mode='medium',seed=seed,quiet=True,verbose=False,dropR2=None)
        e.build()
        node=e.model.nodes['Z']
        if perturb:
            node.Q.params['mean']=node.Q.params['mean']+np.random.default_rng(seed+12307).normal(0,.25,size=node.Q.params['mean'].shape)
            node.updateExpectations()
        initialization_hash=hashlib.sha256(node.getExpectation().tobytes()).hexdigest()
        e.run()
        weights=e.model.nodes['W'].getExpectation();z=e.model.nodes['Z'].getExpectation()
        tau=[np.mean(t,axis=0) if np.ndim(t)==2 else np.asarray(t) for t in e.model.nodes['Tau'].getExpectation()]
        elbo=np.asarray(e.model.getTrainingStats()['elbo']);elbo=elbo[np.isfinite(elbo)]
    return dict(weights=weights,scores=z,precision=tau,elbo=float(elbo[-1]),converged='Converged!' in log.getvalue(),seconds=time.perf_counter()-start,log=log.getvalue(),initialization='PCA plus jitter' if perturb else 'PCA',initialization_hash=initialization_hash)


def tune_mofa(data,seed,rank,config):
    rng=np.random.default_rng(seed+909);masked=[v.copy() for v in data['views']];mask=[]
    for j in [2,3]:
        held=(rng.uniform(size=masked[j].shape)<config['validation_mask_fraction'])&np.isfinite(masked[j]);mask.append(held);masked[j][held]=np.nan
    va,params=fit_preprocess(masked,data['covariates'])
    target=apply_preprocess(data['views'],data['covariates'],params)
    margins=None
    if rank:
        va,margins=rank_transform_fit(va);target=rank_transform_apply(target,margins)
    candidates=[]
    for k in config['factor_grid']:
        fits=[official_mofa(va,k,seed+s,rank,config['mofa_iterations'],perturb=(s==config['starts'][-1])) for s in config['starts']]
        assert len({f['initialization_hash'] for f in fits})==2
        # Starts selected by training evidence, never by masked/test outcomes.
        fit=max(fits,key=lambda f:f['elbo'])
        errors=[]
        for j,held in zip([2,3],mask):
            predicted=fit['scores']@fit['weights'][j].T
            errors.extend((predicted[held]-target[j][held])**2)
        candidates.append(dict(k=k,validation_mse=float(np.mean(errors)),elbo=fit['elbo'],converged=fit['converged'],initializations=[dict(type=f['initialization'],sha256=f['initialization_hash'],elbo=f['elbo'],converged=f['converged']) for f in fits]))
    chosen=min(candidates,key=lambda r:(r['validation_mse'],r['k']))['k']
    full,params=fit_preprocess(data['views'],data['covariates'])
    if rank:full,margins=rank_transform_fit(full)
    fitted=[official_mofa(full,chosen,seed+s,rank,config['mofa_iterations'],perturb=(s==config['starts'][-1])) for s in config['starts']]
    assert len({f['initialization_hash'] for f in fitted})==2
    model=max(fitted,key=lambda f:f['elbo'])
    return dict(model=model,k=chosen,candidates=candidates,params=params,margins=margins,train=full,rank=rank)


def transform(tuned,data):
    views=apply_preprocess(data['views'],data['covariates'],tuned['params'])
    return rank_transform_apply(views,tuned['margins']) if tuned['rank'] else views


def gibbs(x,k,seed,rank,config,multiplier=1):
    n,p=x.shape;missing=~np.isfinite(x);meta=rank_metadata(x);packed=pack_metadata(meta,n,p)
    monitored=[(j,(j+17)%p) for j in range(0,p,3)];correlations=[];traces=[]
    for chain in range(2):
        rng=np.random.default_rng(seed+chain*100003)
        z=np.nan_to_num(x.copy())
        if rank:
            for j,(order,_,_) in enumerate(meta):z[order,j]=ndtri((rankdata(x[order,j])-.5)/len(order))
        h=rng.normal(size=(n,k));l=rng.normal(0,.2,size=(p,k));psi=np.full(p,.5)
        total=np.zeros((p,p));trace=[]
        burn=config['gibbs_burn']*multiplier;draws=config['gibbs_draws']*multiplier
        for iteration in range(burn+draws):
            v=np.linalg.inv(2*np.eye(k)[None,:,:]+(h.T@h)[None,:,:]/psi[:,None,None])
            mu=np.einsum('pk,pkj->pj',z.T@h/psi[:,None],v)
            l=mu+np.einsum('pij,pj->pi',np.linalg.cholesky(v),rng.normal(size=(p,k)))
            v=np.linalg.inv(np.eye(k)+(l.T/psi)@l)
            h=(z/psi)@l@v+rng.normal(size=(n,k))@np.linalg.cholesky(v).T
            psi=(.5+.5*np.sum((z-h@l.T)**2,axis=0))/rng.gamma(2+n/2,size=p)
            if rank:z=latent_update_packed(z,h@l.T,packed,rng,np.sqrt(psi))
            else:
                sample=h@l.T+rng.normal(size=(n,p))*np.sqrt(psi)
                z[missing]=sample[missing]
            if iteration>=burn and (iteration-burn)%config['gibbs_thin']==0:
                c=correlation(l,psi);total+=c;trace.append([c[a,b] for a,b in monitored])
        if rank:assert_order(z,meta)
        correlations.append(total/len(trace));traces.append(trace)
    rh=split_rhat(traces)
    return dict(correlation=np.mean(correlations,axis=0),trace=np.asarray(traces),rhat_max=float(np.max(rh)),extended=multiplier>1)


def fit_gibbs(x,k,seed,rank,config):
    start=time.perf_counter();model=gibbs(x,k,seed,rank,config)
    if model['rhat_max']>=1.1:model=gibbs(x,k,seed,rank,config,multiplier=3)
    model['seconds']=time.perf_counter()-start;model['converged']=model['rhat_max']<1.1
    return model
