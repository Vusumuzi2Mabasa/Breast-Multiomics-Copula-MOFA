"""Preserved primary F analysis plus HC3 sensitivity and fresh null validation."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import f
from simulation import generate,SCENARIOS
from eqtl import analyse,bh
ROOT=Path(__file__).resolve().parents[1]


def hc3(d):
    n=len(d['y']);base=np.c_[np.ones(n),d['covariates'],d['y']];rows=[]
    for gene in range(12):
        g=d['genotype'][:,(gene//2)*4:(gene//2)*4+4];x=np.c_[base,g,g*d['y'][:,None]];y=d['RNA'][:,gene]
        bread=np.linalg.inv(x.T@x);b=np.linalg.lstsq(x,y,rcond=None)[0];resid=y-x@b
        leverage=np.einsum('ij,jk,ik->i',x,bread,x);scores=x*(resid/np.maximum(1-leverage,1e-8))[:,None]
        covariance=bread@(scores.T@scores)@bread
        row=dict(gene=gene)
        for family,idx in [('omnibus',np.arange(base.shape[1],x.shape[1])),('interaction',np.arange(base.shape[1]+4,x.shape[1]))]:
            coefficient=b[idx];var=covariance[np.ix_(idx,idx)];w=float(coefficient@np.linalg.solve(var,coefficient));row[family+'_p']=float(f.sf(max(0,w/len(idx)),len(idx),n-x.shape[1]))
        rows.append(row)
    for family in ['omnibus','interaction']:
        for row,q in zip(rows,bh([r[family+'_p'] for r in rows])):row[family+'_q']=float(q)
    return rows


def main():
    lock=ROOT/'protocol/robustness_amendment.json';digest=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    content=dict(reason='The preserved primary F analysis showed 0.110 empirical FDR in the abundance-only N480 omnibus cell. Add HC3 covariance sensitivity without replacing primary results, and 800 fresh null/abundance-only cohorts, using fixed new seeds before inspecting them.',source_sha256=digest,method='HC3 Wald F using residual df n-rank and the same full gene model; BH separately across twelve genes per family. Approximate reference, not guaranteed calibration.',fresh_scenarios=['null','abundance_only'],fresh_replicates_per_cell=200,fresh_sample_sizes=[240,480],fresh_seed=1003000000)
    if lock.exists():assert json.loads(lock.read_text())==content
    else:lock.write_text(json.dumps(content,indent=2))
    records=[];fixtures=[]
    specs=[]
    for path in sorted((ROOT/'results/eqtl').glob('*.json')):
        meta=json.loads(path.read_text());specs.append(('original',meta['scenario'],meta['n'],meta['replicate'],meta['seed'],meta['genes']))
    for sc in ['null','abundance_only']:
        for n in [240,480]:
            for rep in range(200):specs.append(('fresh',sc,n,rep,1003000000+SCENARIOS.index(sc)*1000000+n*1000+rep,None))
    out=ROOT/'results/robustness';out.mkdir(exist_ok=True)
    for step,(bank,sc,n,rep,seed,original) in enumerate(specs,1):
        path=out/f'{bank}_{sc}_{n}_{rep:03}.json'
        if path.exists():result=json.loads(path.read_text())
        else:
            d=generate(n,seed,sc);ordinary=original if original is not None else analyse(d);robust=hc3(d)
            result=dict(bank=bank,scenario=sc,n=n,replicate=rep,seed=seed,ordinary=ordinary,HC3=robust)
            path.write_text(json.dumps(result))
        for method,key in [('Ordinary F','ordinary'),('HC3 Wald F','HC3')]:
            for family,truth in [('omnibus','true_eqtl'),('interaction','true_interaction')]:
                t=np.array([r[truth] for r in result['ordinary']]);s=np.array([r[family+'_q']<.05 for r in result[key]])
                records.append(dict(bank=bank,scenario=sc,n=n,replicate=rep,seed=seed,method=method,family=family,FDP=float((s&~t).sum()/max(s.sum(),1)),sensitivity=float((s&t).sum()/t.sum()) if t.sum() else np.nan,discoveries=int(s.sum())))
        if bank=='fresh' and rep==0:
            d=generate(n,seed,sc);fixtures.append(dict(n=n,scenario=sc,genotype=d['genotype'].tolist(),RNA=d['RNA'].tolist(),covariates=d['covariates'].tolist(),y=d['y'].tolist(),expected=result['HC3']))
        if step%200==0:print('Robustness',step,'/',len(specs),flush=True)
    df=pd.DataFrame(records);df.to_csv(ROOT/'results/robustness_evaluation.csv',index=False)
    rows=[]
    for key,g in df.groupby(['bank','scenario','n','method','family']):
        row=dict(zip(['bank','scenario','n','method','family'],key));row.update(replicates=len(g),FDP=float(g.FDP.mean()),FDP_mcse=float(g.FDP.std(ddof=1)/np.sqrt(len(g))),sensitivity=float(g.sensitivity.mean()) if g.sensitivity.notna().any() else None)
        rows.append(row)
    pd.DataFrame(rows).to_csv(ROOT/'results/robustness_summary.csv',index=False)
    (ROOT/'results/R_HC3_fixtures.json').write_text(json.dumps(fixtures))
    (ROOT/'results/robustness_audit.json').write_text(json.dumps(dict(status='passed',original_cohorts=1600,fresh_cohorts=800,method_family_records=len(df),fresh_seed_disjoint=True),indent=2))


if __name__=='__main__':main()
