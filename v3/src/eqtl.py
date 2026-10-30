"""Gene-level cis association and interaction tests with independent truth."""
import numpy as np
from scipy.stats import f


def bh(p):
    p=np.asarray(p);order=np.argsort(p);ranked=p[order]*len(p)/np.arange(1,len(p)+1)
    result=np.empty(len(p));result[order]=np.minimum(1,np.minimum.accumulate(ranked[::-1])[::-1]);return result


def fit(x,y):
    coefficient,_,rank,_=np.linalg.lstsq(x,y,rcond=None)
    residual=y-x@coefficient
    return float(residual@residual),int(rank),coefficient


def nested(x0,x1,y):
    a,r0,_=fit(x0,y);b,r1,coef=fit(x1,y)
    df=r1-r0;den=len(y)-r1
    if df<=0 or den<=0:raise ValueError('Invalid model rank')
    statistic=max(0,(a-b)/df)/(b/den)
    return float(f.sf(statistic,df,den)),statistic,coef


def analyse(d):
    n=len(d['y']);base=np.c_[np.ones(n),d['covariates'],d['y']]
    rows=[]
    for j in range(12):
        indices=np.arange((j//2)*4,(j//2)*4+4);g=d['genotype'][:,indices]
        common=np.c_[base,g];full=np.c_[common,g*d['y'][:,None]];y=d['RNA'][:,j]
        p,stat,coef=nested(base,full,y);pi,stati,_=nested(common,full,y)
        tags=[]
        for snp in indices:
            ps,_,cs=nested(base,np.c_[base,d['genotype'][:,snp]],y)
            tags.append((ps,int(snp),float(cs[-1])))
        lead=min(tags)
        causal=(0 if j%2==0 else 2)
        rows.append(dict(gene=j,omnibus_p=p,omnibus_F=stat,interaction_p=pi,interaction_F=stati,causal_main_estimate=float(coef[base.shape[1]+causal]),causal_interaction_estimate=float(coef[base.shape[1]+4+causal]),lead_variant=lead[1],lead_p=lead[0],lead_effect=lead[2],true_eqtl=bool(d['beta'][j]!=0 or d['delta'][j]!=0),true_interaction=bool(d['delta'][j]!=0)))
    for field in ['omnibus','interaction']:
        for row,q in zip(rows,bh([r[field+'_p'] for r in rows])):row[field+'_q']=float(q)
    return rows
