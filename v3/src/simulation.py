"""Small, wholly synthetic genotype/RNA breast cancer case cohorts."""
from pathlib import Path
import json
import numpy as np
from scipy.special import expit, ndtr

ROOT=Path(__file__).resolve().parents[1]
SCENARIOS=['null','shared_eqtl','subtype_eqtl','abundance_only']


def generate(n,seed,scenario):
    if scenario not in SCENARIOS:raise ValueError(scenario)
    rng=np.random.default_rng(seed)
    ancestry=rng.choice(3,n,p=[.5,.3,.2])
    age=rng.normal(0,1,n);cells=rng.beta(5,5,n);batch=rng.integers(0,3,n)
    cov=np.c_[age,(cells-.5)/.15,ancestry==1,ancestry==2,batch==1,batch==2].astype(float)
    g=np.empty((n,24));ld=.65**np.abs(np.arange(4)[:,None]-np.arange(4)[None,:])
    maf=np.tile([.18,.25,.33,.4],6)
    for locus in range(6):
        u=ndtr(rng.multivariate_normal(np.zeros(4),ld,size=n));f=np.clip(maf[locus*4:locus*4+4]+(ancestry[:,None]-0.7)*.025,.08,.46)
        g[:,locus*4:locus*4+4]=(u>(1-f)**2).astype(int)+(u>1-f*f).astype(int)
    y=rng.binomial(1,expit(-.65+.3*age+.25*cov[:,1]+.2*(ancestry==2)+.25*(g[:,0]-2*maf[0])))
    h=rng.normal(size=(n,3));loading=np.zeros((12,3))
    for gene in range(12):loading[gene,gene//4]=.25
    beta=np.zeros(12);delta=np.zeros(12);gamma=np.zeros(12)
    if scenario in ['shared_eqtl','subtype_eqtl']:beta[::2]=.45
    if scenario=='subtype_eqtl':delta[[0,4,8]]=.65;delta[[3,7,11]]=-.65
    if scenario in ['subtype_eqtl','abundance_only']:gamma[[0,1,4,5,8,9]]=.65
    latent=np.empty((n,12))
    for gene in range(12):
        causal=(gene//2)*4+(0 if gene%2==0 else 2)
        latent[:,gene]=5+.2*age+.2*cov[:,1]+.12*cov[:,2]+.16*cov[:,3]+.15*cov[:,4]-.15*cov[:,5]+h@loading[gene]+rng.normal(0,.35,n)+g[:,causal]*beta[gene]+g[:,causal]*y*delta[gene]+y*gamma[gene]
    library=rng.lognormal(np.log(20e6),.2,n)
    mu=library[:,None]/1e6*2**latent
    counts=rng.negative_binomial(30,30/(30+mu))
    rna=np.log2(counts/library[:,None]*1e6+1)
    return dict(genotype=g,RNA=rna,counts=counts,covariates=cov,y=y,ancestry=ancestry,batch=batch,age=50+10*age,cells=cells,library=library,h=h,beta=beta,delta=delta,gamma=gamma,seed=seed,scenario=scenario)


def export_example(path=None):
    import pandas as pd
    out=path or ROOT/'data/example';out.mkdir(parents=True,exist_ok=True)
    d=generate(240,805240000,'subtype_eqtl')
    for name in ['genotype','RNA','counts','covariates']:
        cols=[f'SIM_SNP_{j+1:02}' for j in range(24)] if name=='genotype' else ([f'SIM_GENE_{j+1:02}' for j in range(12)] if name in ['RNA','counts'] else ['age_z','cell_fraction_z','ancestry_B','ancestry_C','batch_B','batch_C'])
        pd.DataFrame(d[name],columns=cols).to_csv(out/f'{name}.csv',index=False)
    pd.DataFrame(dict(id=[f'SIM_{j+1:03}' for j in range(240)],TNBC=d['y'],ancestry=d['ancestry'],batch=d['batch'],age=d['age'],cell_fraction=d['cells'],library_size=d['library'])).to_csv(out/'metadata.csv',index=False)
    genes=[];variants=[]
    for locus in range(6):
        for j in range(4):variants.append(dict(variant=f'SIM_SNP_{4*locus+j+1:02}',synthetic_chromosome=locus+1,position=100000+j*10000,maf_target=float(np.tile([.18,.25,.33,.4],6)[4*locus+j])))
        for j in range(2):
            gene=2*locus+j;genes.append(dict(gene=f'SIM_GENE_{gene+1:02}',synthetic_chromosome=locus+1,TSS=115000+j*10000,causal_variant=f'SIM_SNP_{4*locus+(0 if j==0 else 2)+1:02}',main_beta=float(d['beta'][gene]),interaction_beta=float(d['delta'][gene]),subtype_abundance_shift=float(d['gamma'][gene])))
    pd.DataFrame(genes).to_csv(out/'gene_truth.csv',index=False);pd.DataFrame(variants).to_csv(out/'variant_map.csv',index=False)
    (out/'manifest.json').write_text(json.dumps(dict(disclosure='Entirely synthetic breast cancer cases, genomic positions and gene IDs. Not real breast cancer susceptibility or eQTL findings.',seed=d['seed'],scenario=d['scenario'],samples=240,views=['genotype','RNA'],RNA_units='log2 CPM plus one using known simulated whole-assay library',genotype_units='0/1/2 alternate allele dosage',case_only=True),indent=2))


if __name__=='__main__':export_example()
