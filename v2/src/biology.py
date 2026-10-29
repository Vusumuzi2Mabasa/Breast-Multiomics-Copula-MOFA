"""Synthetic gene modules with distinct genomic, RNA and protein measurements.

No real ancestry reference, participant, gene annotation or tumour data are used.
The generating order is a simulation mechanism, not causal inference from data.
"""
import numpy as np
from scipy.special import expit,ndtr
from scipy.stats import binom,nbinom,t

NULLS=['null','abundance_only','scale_only','technical_null','missing_null']
ALTERNATIVES=['coordination_moderate','coordination_strong','both','opposing_coordination','missing_coordination']
STRESS=['heavy_tail_null','heavy_tail_coordination']
GENES=16


def generate(n,seed,scenario):
    rng=np.random.default_rng(seed)
    country=rng.choice(5,n,p=[.35,.25,.20,.12,.08])
    age=rng.normal(size=n);cells=rng.normal(size=n)
    treatment=rng.binomial(1,expit(-.8+.25*age+.15*country))
    batch=(country+rng.integers(0,3,n))%3
    # Four two-variant LD blocks, entirely synthetic; no African-LD claim.
    block=rng.normal(size=(n,4))
    zg=.7*np.repeat(block,2,axis=1)+np.sqrt(1-.7**2)*rng.normal(size=(n,8))
    maf=np.clip(np.linspace(.12,.36,8)[None,:]+.015*(country[:,None]-1.4),.05,.45)
    g=binom.ppf(np.clip(ndtr(zg),1e-8,1-1e-8),2,maf)
    y=rng.binomial(1,expit(-1+.25*age+.45*cells+.18*(country-1.4)+.2*(g[:,0]-2*maf[:,0])))
    somprob=np.linspace(-3,-1.8,8)[None,:]+.15*age[:,None]
    if scenario in ['abundance_only','both','technical_null','missing_null','missing_coordination']:
        somprob=np.broadcast_to(somprob,(n,8)).copy();somprob[:,:4]+=.4*y[:,None]
    som=rng.binomial(1,expit(somprob))
    h=rng.normal(size=(n,4))
    zr=.45*np.repeat(h,4,axis=1)+np.sqrt(1-.45**2)*rng.normal(size=(n,16))
    eta0=np.arctanh(.2)+.25*cells+.15*age+.08*(country-1.4)
    beta=np.zeros(16)
    if scenario in ['coordination_moderate','both','missing_coordination']:beta[:4]=.35
    if scenario in ['coordination_strong','heavy_tail_coordination']:beta[:4]=.65
    if scenario=='opposing_coordination':beta[:4]=[.45,-.45,.45,-.45]
    rho=np.tanh(eta0[:,None]+y[:,None]*beta)
    zp=rho*zr+np.sqrt(1-rho*rho)*rng.normal(size=(n,16))
    if scenario.startswith('heavy_tail'):
        scale=np.sqrt(rng.chisquare(5,(n,1))/5);zr/=scale;zp/=scale
    mean_r=np.zeros((n,16));mean_p=np.zeros((n,16))
    for j in range(16):
        cis=j%8
        mean_r[:,j]=.30*(g[:,cis]-2*maf[:,cis])+.25*som[:,cis]+.15*age+.25*cells+.12*treatment
        mean_p[:,j]=(.22 if j%4!=3 else -.15)*(g[:,cis]-2*maf[:,cis])+.20*som[:,cis]+.12*age+.30*cells+.10*treatment
    shifted=scenario in ['abundance_only','both','technical_null','missing_null','missing_coordination']
    if shifted:
        mean_r[:,:4]+=.45*y[:,None];mean_p[:,:4]+=.30*y[:,None]
    if scenario in ['technical_null','missing_null','missing_coordination']:
        mean_r+=.18*country[:,None]+.18*batch[:,None]
        mean_p+=.12*country[:,None]+.25*batch[:,None]
    sr=np.ones((n,1));sp=np.ones((n,1))
    if scenario in ['scale_only','technical_null','missing_null','missing_coordination']:
        sr=np.exp(.20*y[:,None]+.12*cells[:,None]);sp=np.exp(-.15*y[:,None]+.12*age[:,None])
    # Conditional NB margins, with a known full-assay library size, followed by
    # routine log2 CPM; counts are not obtained by summing the retained genes.
    library=rng.lognormal(np.log(2e7),.25,n)
    mu=np.exp(4.5+mean_r)*library[:,None]/2e7
    u=t.cdf(zr,5) if scenario.startswith('heavy_tail') else ndtr(zr)
    # Scale effects operate on the latent RNA quantile; regular variance = 1.
    if not scenario.startswith('heavy_tail'):u=ndtr(zr*sr)
    counts=nbinom.ppf(np.clip(u,1e-8,1-1e-8),12,12/(12+mu))
    rna=np.log2(counts/library[:,None]*1e6+1)
    protein_log=3+mean_p+sp*zp
    missing=np.zeros_like(protein_log,dtype=bool)
    if scenario in ['missing_null','missing_coordination']:
        prob=expit(-2.4-.65*(protein_log-np.mean(protein_log,axis=0)))
        missing=rng.uniform(size=protein_log.shape)<prob
    protein=protein_log.copy();protein[missing]=np.nan
    # Common nuisance terms; genotype and somatic adjustment is pair-specific.
    x=np.c_[age,cells,treatment,np.eye(5)[country][:,1:],np.eye(3)[batch][:,1:]]
    return dict(germline=g,somatic=som,RNA=rna,protein=protein,
                views=[g,som,rna,protein],y=y,country=country,batch=batch,
                covariates=x,h=h,latent_RNA=zr,latent_protein=zp,
                rho_null=np.tanh(eta0),beta=beta,active=np.arange(16)<4 if np.any(beta) else np.zeros(16,bool),
                counts=counts,library=library,missing=missing,seed=seed,scenario=scenario)


def pair_design(d,j):
    x=d['covariates'];y=d['y'];cis=j%8
    mean=np.c_[np.ones(len(y)),y,x,d['germline'][:,cis],d['somatic'][:,cis],y*x[:,0],y*x[:,1]]
    scale=np.c_[np.ones(len(y)),y,x[:,:2]]
    dependence=np.c_[np.ones(len(y)),x]
    return mean,scale,dependence
