"""Conditional coordination statistics; family calibration occurs separately."""
import time
import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata
from scipy.optimize import minimize
from biology import pair_design

METHODS=['Partial Pearson','Partial Spearman','Normal-score difference','Product regression','Conditional copula','Oracle margins']


def residualise(y,x):
    b=np.linalg.lstsq(x,y,rcond=None)[0]
    return y-x@b


def crossfit_margin(y,mean,scale,fold):
    out=np.full(len(y),np.nan)
    for k in [0,1]:
        train=(fold!=k)&np.isfinite(y);test=(fold==k)&np.isfinite(y)
        b=np.linalg.lstsq(mean[train],y[train],rcond=None)[0]
        e=y[train]-mean[train]@b
        floor=.1*np.mean(e*e)+1e-6
        a=np.linalg.lstsq(scale[train],np.log(e*e+floor),rcond=None)[0]
        logvar=np.clip(scale[test]@a,-6,6)
        out[test]=(y[test]-mean[test]@b)/np.exp(.5*logvar)
    return out


def score_normalise(x,t):
    z=np.zeros_like(x)
    for g in [0,1]:
        take=t==g;n=sum(take)
        z[take]=ndtri((rankdata(x[take])-.5)/n)
    return z


def correlation_stat(a,b,t):
    r=[];sizes=[]
    for group in [0,1]:
        take=t==group;sizes.append(sum(take))
        r.append(np.clip(np.corrcoef(a[take],b[take])[0,1],-.999,.999))
    return float((np.arctanh(r[1])-np.arctanh(r[0]))/np.sqrt(1/(sizes[0]-3)+1/(sizes[1]-3)))


def copula_objective(theta,u,v,x):
    eta=x@theta;rho=np.tanh(eta);one=np.maximum(1-rho*rho,1e-8)
    a=u*u+v*v;b=u*v
    ll=-.5*np.log(one)-.5*(a-2*rho*b)/one
    s=rho+((1+rho*rho)*b-rho*a)/one
    return -float(np.mean(ll)),-x.T@s/len(u)


def copula_score(u,v,t,x):
    start=np.zeros(x.shape[1]);start[0]=np.arctanh(np.clip(np.corrcoef(u,v)[0,1],-.8,.8))
    fit=minimize(copula_objective,start,args=(u,v,x),jac=True,method='L-BFGS-B',bounds=[(-2,2)]*len(start),
                 options={'maxiter':150,'ftol':1e-9,'gtol':1e-5})
    value,gradient=copula_objective(fit.x,u,v,x)
    if not fit.success and np.max(np.abs(gradient))>1e-3:raise RuntimeError('Conditional copula null fit did not converge')
    rho=np.tanh(x@fit.x);one=np.maximum(1-rho*rho,1e-8)
    s=rho+((1+rho*rho)*u*v-rho*(u*u+v*v))/one
    w=1+rho*rho
    b=np.linalg.lstsq(x*np.sqrt(w[:,None]),t*np.sqrt(w),rcond=None)[0]
    efficient=t-x@b
    contrib=efficient*s
    stat=np.sum(contrib)/np.sqrt(np.sum(contrib**2)+1e-12)
    return float(stat),float(np.max(np.abs(gradient)))


def product_score(u,v,t,x):
    design=np.c_[x,t,x[:,1]**2,x[:,2]**2]
    inv=np.linalg.pinv(design.T@design)
    b=inv@design.T@(u*v);e=u*v-design@b
    leverage=np.sum((design@inv)*design,axis=1)
    meat=design.T@((e/np.maximum(1-leverage,.05))[:,None]**2*design)
    var=inv@meat@inv
    col=x.shape[1]
    return float(b[col]/np.sqrt(max(var[col,col],1e-12)))


def evaluate(d):
    n=len(d['y']);rng=np.random.default_rng(d['seed']+70000000)
    folds=np.arange(n)%2;rng.shuffle(folds)
    output=np.zeros((len(METHODS),16));gradient=np.zeros(16);sizes=[]
    start=time.perf_counter()
    for j in range(16):
        mean,scale,x=pair_design(d,j)
        r=d['RNA'][:,j];p=d['protein'][:,j]
        good=np.isfinite(r)&np.isfinite(p);t=d['y'][good]
        sizes.append(int(sum(good)))
        if min(sum(t==0),sum(t==1))<12:raise RuntimeError('Insufficient subtype observations')
        a=residualise(r[good],mean[good]);b=residualise(p[good],mean[good])
        output[0,j]=correlation_stat(a,b,t)
        ar=residualise(rankdata(r[good]),mean[good]);br=residualise(rankdata(p[good]),mean[good])
        output[1,j]=correlation_stat(ar,br,t)
        a=crossfit_margin(r,mean,scale,folds)[good];b=crossfit_margin(p,mean,scale,folds)[good]
        u=score_normalise(a,t);v=score_normalise(b,t)
        output[2,j]=correlation_stat(u,v,t)
        output[3,j]=product_score(u,v,t,x[good])
        output[4,j],gradient[j]=copula_score(u,v,t,x[good])
        # Oracle has access to latent pre-measurement scores unavailable to a
        # real analysis. Missing pairs use the same selected participants.
        output[5,j],_=copula_score(d['latent_RNA'][good,j],d['latent_protein'][good,j],t,x[good])
    return dict(statistics=output,gradient=gradient,pair_sizes=np.array(sizes),seconds=time.perf_counter()-start)
