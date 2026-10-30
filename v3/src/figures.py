"""Readable simulation figures, example QC and data-driven paper tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from simulation import generate
from eqtl import analyse,nested,fit
ROOT=Path(__file__).resolve().parents[1]
SC=['null','shared_eqtl','subtype_eqtl','abundance_only']
SL=['Null','Shared eQTL','Subtype eQTL','Abundance only']
METHOD=['MOFA+','DIABLO','PCA','Combined ridge','Genotype ridge','RNA ridge','Clinical only']
COLOR=['#2166ac','#b2182b','#5e3c99','#1b7837','#e08214','#4393c3','#777777']


def save(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    p=ROOT/f'figures/{name}.svg';fig.savefig(p,bbox_inches='tight',facecolor='white');plt.close(fig)
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n',encoding='utf-8',newline='\n')


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.right':False,'axes.spines.top':False})
    i=pd.read_csv(ROOT/'results/integration_summary.csv',keep_default_na=False,na_values=['']);q=pd.read_csv(ROOT/'results/eqtl_summary.csv',keep_default_na=False,na_values=[''])
    example=generate(240,805240000,'subtype_eqtl');rows=analyse(example)
    table=pd.DataFrame(rows);table.to_csv(ROOT/'results/example_eqtl.csv',index=False)
    # Genotype QC is descriptive; HWE is evaluated within synthetic ancestry strata.
    g=example['genotype'];freq=g.mean(0)/2;maf=np.minimum(freq,1-freq);hwe=[]
    for anc in range(3):
        for j in range(24):
            x=g[example['ancestry']==anc,j];p=x.mean()/2;expected=len(x)*np.array([(1-p)**2,2*p*(1-p),p*p]);obs=np.bincount(x.astype(int),minlength=3)
            reliable=bool(np.all(expected>=5));pv=float(chi2.sf(np.sum((obs-expected)**2/np.maximum(expected,1e-12)),1)) if reliable else None
            hwe.append(dict(ancestry=anc,variant=j,samples=len(x),minimum_expected=float(expected.min()),approximation_reliable=reliable,p=pv))
    pd.DataFrame(hwe).to_csv(ROOT/'results/example_HWE.csv',index=False)
    qc=dict(status='passed',samples=240,genotype_variants=24,RNA_genes=12,genotype_call_rate=float(np.isfinite(g).mean()),minimum_MAF=float(maf.min()),genotype_values_valid=bool(np.isin(g,[0,1,2]).all()),finite_RNA=bool(np.isfinite(example['RNA']).all()),nonnegative_integer_counts=bool(np.all(example['counts']>=0) and np.all(example['counts']==np.floor(example['counts']))),positive_libraries=bool(np.all(example['library']>0)),TNBC=int(example['y'].sum()),non_TNBC=int((1-example['y']).sum()),protein_features=0,somatic_features=0,synthetic_ancestry_counts=np.bincount(example['ancestry'],minlength=3).tolist(),hwe_reliable_cells=sum(x['approximation_reliable'] for x in hwe),hwe_policy='Pearson approximation only where all expected genotype counts >=5, within ancestry; descriptive, no phenotype-dependent filtering; no real genotype QC validation implied.')
    assert qc['genotype_values_valid'] and qc['finite_RNA'] and qc['nonnegative_integer_counts'] and qc['positive_libraries'] and qc['minimum_MAF']>.05
    (ROOT/'results/example_QC.json').write_text(json.dumps(qc,indent=2))
    fig,ax=plt.subplots(figsize=(7.2,4.0));ax.axis('off');ax.set_xlim(0,10);ax.set_ylim(0,5)
    ax.text(2,4.3,'Germline genotype\n24 variants in 6 synthetic loci',ha='center',fontsize=12)
    ax.text(7.6,4.3,'RNA expression\n12 synthetic genes',ha='center',fontsize=12)
    ax.annotate('cis eQTL and subtype interaction tests',xy=(5,3.1),xytext=(5,3.1),ha='center',fontsize=11,color='#2166ac')
    for x in [2,7.6]:ax.annotate('',xy=(5,3.5),xytext=(x,4),arrowprops=dict(arrowstyle='->',color='#777777'))
    ax.text(2.6,2,'MOFA+\nUnsupervised factors',ha='center',fontsize=12)
    ax.text(7.5,2,'DIABLO\nSupervised subtype signature',ha='center',fontsize=12)
    ax.text(5,.7,'Independent TNBC versus non TNBC prediction\nClinical and technical covariates available to both workflows',ha='center',fontsize=10)
    for x in [2.6,7.5]:ax.annotate('',xy=(5,1.35),xytext=(x,1.6),arrowprops=dict(arrowstyle='->',color='#777777'))
    save(fig,'study_design')
    fig,axes=plt.subplots(2,1,figsize=(7.2,6.2),sharex=True)
    for ax,n in zip(axes,[240,480]):
        for j,(m,c) in enumerate(zip(METHOD,COLOR)):
            r=i[(i.n==n)&(i.method==m)].set_index('scenario').loc[SC]
            ax.errorbar(np.arange(4)+(j-3)*.075,r.balanced_accuracy,yerr=1.96*r.balanced_accuracy_mcse,fmt='o',ms=3,capsize=1.5,color=c,label=m)
        ax.axhline(.5,color='grey',ls=':');ax.set(title=f'Independent test balanced accuracy, training N={n}',ylabel='Balanced accuracy',ylim=(.4,1.02),xticks=np.arange(4),xticklabels=SL)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=4,bbox_to_anchor=(.5,1.06),fontsize=8.5)
    fig.tight_layout();save(fig,'integration_performance')
    fig,axes=plt.subplots(2,2,figsize=(7.2,6.1),sharex=True)
    for col,n in enumerate([240,480]):
        for j,(family,label,color) in enumerate([('omnibus','Any cis association','#2166ac'),('interaction','Subtype interaction','#b2182b')]):
            r=q[(q.n==n)&(q.family==family)].set_index('scenario').loc[SC]
            for row,(metric,ylabel) in enumerate([('sensitivity','Sensitivity among true genes'),('FDP','Mean false discovery proportion')]):
                mu=r[metric].to_numpy();err=1.96*r[metric+'_mcse'].to_numpy()
                axes[row,col].errorbar(np.arange(4)+(j-.5)*.12,mu,yerr=err,fmt='o',capsize=2,color=color,label=label)
                axes[row,col].set(xticks=np.arange(4),xticklabels=SL,ylabel=ylabel,title=f'N={n}',ylim=(-.025,1.025) if row==0 else (-.01,max(.12,float(q.FDP.max()+2*q.FDP_mcse.max()))))
        axes[1,col].axhline(.05,color='grey',ls=':')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=2,bbox_to_anchor=(.5,1.04),fontsize=9)
    for ax in axes[1]:ax.tick_params(axis='x',rotation=25)
    fig.tight_layout();save(fig,'eqtl_performance')
    robust=pd.read_csv(ROOT/'results/robustness_summary.csv',keep_default_na=False,na_values=[''])
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.7),sharey=True)
    for ax,n in zip(axes,[240,480]):
        for j,(bank,method,color,marker) in enumerate([('original','Ordinary F','#2166ac','o'),('original','HC3 Wald F','#b2182b','o'),('fresh','Ordinary F','#2166ac','s'),('fresh','HC3 Wald F','#b2182b','s')]):
            r=robust[(robust.n==n)&(robust.bank==bank)&(robust.method==method)&(robust.family=='omnibus')].set_index('scenario').loc[['null','abundance_only']]
            ax.errorbar(np.arange(2)+(j-1.5)*.1,r.FDP,yerr=1.96*r.FDP_mcse,fmt=marker,color=color,ms=4,capsize=2,label=f'{method}, {bank}')
        ax.axhline(.05,color='grey',ls=':');ax.set(title=f'No genetic expression effects, N={n}',xticks=[0,1],xticklabels=['Null','Abundance only'],ylabel='Any false cis gene',ylim=(-.01,max(.18,float(robust.FDP.max()+2*robust.FDP_mcse.max()))))
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.09),fontsize=8.5)
    fig.tight_layout();save(fig,'eqtl_robustness')
    fixture=np.load(ROOT/'results/integration/subtype_eqtl_240_00.npz')
    diab=json.loads((ROOT/'results/diablo/subtype_eqtl_240_00.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.3));z=fixture['scores'];y=example['y'];dv=diab['train_scores'];zg=np.asarray(dv['genotype']);zr=np.asarray(dv['RNA']);zg=zg.reshape(len(y),-1);zr=zr.reshape(len(y),-1)
    for group,color,label in [(0,'#2166ac','Non TNBC'),(1,'#b2182b','TNBC')]:
        idx=y==group;axes[0].scatter(z[idx,0],z[idx,1],s=12,alpha=.6,color=color,label=label);axes[1].scatter(zg[idx,0],zr[idx,0],s=12,alpha=.6,color=color)
    axes[0].set(title='MOFA+ training representation',xlabel='Factor 1',ylabel='Factor 2');axes[1].set(title='DIABLO training representation',xlabel='Genotype component 1',ylabel='RNA component 1');axes[0].legend(fontsize=8)
    fig.tight_layout();save(fig,'example_factors')
    base=np.c_[np.ones(240),example['covariates'],y];gg=example['genotype'][:,:4];full=np.c_[base,gg,gg*y[:,None]]
    _,_,coef=fit(full,example['RNA'][:,0]);adjusted=example['RNA'][:,0]-example['covariates']@coef[1:7]
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.5));rng=np.random.default_rng(951)
    for group,color,label in [(0,'#2166ac','Non TNBC'),(1,'#b2182b','TNBC')]:
        idx=y==group;axes[0].scatter(gg[idx,0]+rng.normal(0,.035,idx.sum()),adjusted[idx],color=color,s=10,alpha=.5,label=label)
        means=[np.mean(adjusted[(y==group)&(gg[:,0]==v)]) if np.any((y==group)&(gg[:,0]==v)) else np.nan for v in range(3)]
        axes[0].plot(range(3),means,'o-',color=color)
    axes[0].set(title='Illustrative genotype expression pattern',xlabel='SIM SNP 01 alternate allele dosage',ylabel='Covariate adjusted log2 CPM',xticks=[0,1,2]);axes[0].legend(fontsize=8)
    # Marginal cis tag associations: distinct from gene-level joint/interaction tests.
    ps=[nested(base,np.c_[base,gg[:,j]],example['RNA'][:,0])[0] for j in range(4)]
    axes[1].plot([100,110,120,130],-np.log10(np.maximum(ps,1e-300)),'o-',color='#2166ac');axes[1].set(title='SIM GENE 01 local association',xlabel='Synthetic position in kb',ylabel='Negative log10 nominal P')
    fig.tight_layout();save(fig,'example_eqtl')
    fig,ax=plt.subplots(figsize=(5.2,4.7));matrix=-np.log10(np.maximum(table[['omnibus_q','interaction_q']].to_numpy(),1e-15));im=ax.imshow(matrix,aspect='auto',cmap='Blues',vmin=0,vmax=max(3,float(matrix.max())))
    ax.set(xticks=[0,1],xticklabels=['Any cis association','Subtype interaction'],yticks=np.arange(12),yticklabels=[f'SIM GENE {j+1:02}' for j in range(12)])
    for row in range(12):
        for col,name in enumerate(['omnibus_q','interaction_q']):
            if table.iloc[row][name]<.05:ax.text(col,row,'*',ha='center',va='center',color='white' if matrix[row,col]>matrix.max()*.6 else 'black',fontsize=14)
    fig.colorbar(im,ax=ax,label='Negative log10 BH adjusted P');fig.tight_layout();save(fig,'example_qtl_map')
    print('Seven figures and independent example QC generated')


if __name__=='__main__':main()
