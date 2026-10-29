"""Small synthetic multi-omics panels; no participant or collaborator data."""
import numpy as np
from scipy.special import ndtr, expit
from scipy.stats import binom, nbinom, t

SIZES=[6,6,24,12]
NAMES=['germline','somatic','RNA','protein']
SCENARIOS=['regular','skewed','technical','missing','heavy_tail','dependence_only']


def generate(n, seed, scenario):
    rng=np.random.default_rng(seed)
    country=rng.choice(5,n,p=[.35,.25,.20,.12,.08])
    age=rng.normal(size=n)
    cells=rng.normal(size=n)
    batch=rng.integers(0,3,n)
    h=rng.normal(size=(n,3))
    p=sum(SIZES)
    # Two shared biological programmes and a weaker shared nuisance programme.
    load=np.zeros((p,3))
    for start,size in zip(np.cumsum([0]+SIZES[:-1]),SIZES):
        for local in range(size):
            load[start+local,local%2]=(.9 if local<min(size,8) else .45)*(-1 if local%5==0 else 1)
            load[start+local,2]=.25
    eps=rng.normal(size=(n,p))
    latent=h@load.T+eps
    sd=np.sqrt(1+np.sum(load**2,axis=1))
    latent/=sd
    truth=(load@load.T+np.eye(p))/sd[:,None]/sd[None,:]
    if scenario=='heavy_tail':
        latent/=np.sqrt(rng.chisquare(4,n)[:,None]/4)
        u=t.cdf(latent,df=4)
    else:u=ndtr(latent)
    y=rng.binomial(1,expit(-1+.9*h[:,0]+.5*h[:,1]+.15*age+.12*(country-1.4)))
    if scenario=='dependence_only':
        y=rng.binomial(1,.30,n)
        # Same standard normal feature margins in both subtypes; only a
        # prespecified RNA/protein pair changes correlation (0 versus .65).
        rho=.65*y
        latent[:,36]=rho*latent[:,12]+np.sqrt(1-rho*rho)*rng.normal(size=n)
        u[:,36]=ndtr(latent[:,36])
        # No subtype mean shifts anywhere. Pair test is separate from factors.
    u=np.clip(u,1e-8,1-1e-8)
    maf=np.linspace(.12,.4,6)
    germline=binom.ppf(u[:,:6],2,maf)
    somatic=(u[:,6:12]>1-np.linspace(.04,.20,6)).astype(float)
    library=rng.lognormal(np.log(2e7),.35,n)
    if scenario=='regular':
        rna=7+latent[:,12:36]
        protein=10+latent[:,36:48]
        raw_counts=None
    else:
        mean=np.exp(np.linspace(3,6,24))[None,:]*(library[:,None]/2e7)
        counts=nbinom.ppf(u[:,12:36],8,8/(8+mean))
        raw_counts=counts.copy()
        rna=np.log2(counts/library[:,None]*1e6+1)
        abundance=np.exp(2+.7*latent[:,36:48])
        protein=np.log2(abundance)
    if scenario in ['technical','missing']:
        # Overlapping, imperfect confounding; no correction can infer away
        # completely aliased country and biological effects.
        rna+=.35*cells[:,None]+.2*(country[:,None]-1.4)+.15*batch[:,None]
        protein+=.45*cells[:,None]+.15*(country[:,None]-1.4)+.2*batch[:,None]
    if scenario=='missing':
        pr=expit(-2-.8*(protein-np.mean(protein,axis=0)))
        protein[rng.uniform(size=protein.shape)<pr]=np.nan
    cov=np.column_stack([age, np.eye(5)[country][:,1:], cells, np.eye(3)[batch][:,1:]])
    return dict(views=[germline,somatic,rna,protein],y=y,h=h,
                covariates=cov,country=country,truth=truth,raw_counts=raw_counts,
                seed=seed,scenario=scenario,loadings=load)


def standardise(train,test):
    """Train-only scaling and missing-value means; binary view untransformed."""
    a=[];b=[]
    for i,(x,z) in enumerate(zip(train,test)):
        if i==1:a.append(x.copy());b.append(z.copy());continue
        mu=np.nanmean(x,axis=0);sd=np.nanstd(x,axis=0)
        sd=np.maximum(sd,.1)
        a.append((x-mu)/sd);b.append((z-mu)/sd)
    return a,b
