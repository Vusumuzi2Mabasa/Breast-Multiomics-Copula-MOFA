"""Summarise all fixed evaluation replicates and create publication figures."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm
from simulate import SCENARIOS,generate
from models import normal_scores_fit
from metrics import full_subspace_recovery,basis
ROOT=Path(__file__).resolve().parents[1]
METHODS=['MOFA+','MOFA+ Gaussian','Copula','Rank PCA','Clinical only','Oracle factors']
COLORS={'MOFA+':'#2166ac','MOFA+ Gaussian':'#67a9cf','Copula':'#b2182b','Rank PCA':'#5e3c99','Clinical only':'#999999','Oracle factors':'#1b7837'}
LABELS={'regular':'Regular margins','skewed':'Skewed measurements','technical':'Technical variation','missing':'Missing proteins','heavy_tail':'Heavy-tailed dependence','dependence_only':'Dependence change only'}


def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half


def savefig(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    fig.savefig(ROOT/f'figures/{name}.svg',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    protocol=json.loads((ROOT/'protocol/design.json').read_text())
    paths=sorted((ROOT/'results/replicates').glob('*.json'))
    assert len(paths)==len(SCENARIOS)*2*protocol['replicates'],len(paths)
    rows=[r for p in paths for r in json.loads(p.read_text())]
    df=pd.DataFrame(rows)
    df['original_rank_conditional_recovery']=df['recovery']
    df['estimated_rank']=np.nan
    df['predictive_recovery']=np.nan
    for (scenario,n,rep),g in df.groupby(['scenario','n','replicate']):
        f=np.load(ROOT/f'results/replicates/{scenario}_{n}_{rep}_fixture.npz')
        ns,_=normal_scores_fit(f['data']);_,_,v=np.linalg.svd(ns,full_matrices=False)
        rank_scores=ns@v[:3].T
        for idx,row in g.iterrows():
            method=row['method']
            if method=='Clinical only' or row['status']!='ok':continue
            key={'MOFA+':'mofa_scores','MOFA+ Gaussian':'mofa_gaussian_scores','Copula':'copula_scores','Oracle factors':'h'}.get(method)
            scores=rank_scores if method=='Rank PCA' else f[key]
            df.loc[idx,'recovery']=full_subspace_recovery(scores,f['h'])
            df.loc[idx,'estimated_rank']=basis(scores).shape[1]
            df.loc[idx,'predictive_recovery']=full_subspace_recovery(scores,f['h'][:,:2])
    df.to_csv(ROOT/'results/evaluation.csv',index=False)
    df[['scenario','n','replicate','method','original_rank_conditional_recovery','estimated_rank','recovery']].to_csv(ROOT/'results/recovery_correction.csv',index=False)
    assert len(df)==len(paths)*len(METHODS)
    assert not df.duplicated(['scenario','n','replicate','method']).any()
    ok=df[df.status=='ok'];summary=[]
    for (s,n,m),g in ok.groupby(['scenario','n','method'],sort=False):
        item=dict(scenario=s,n=n,method=m,successful=len(g),planned=protocol['replicates'])
        for metric in ['auc','ap','brier','recovery','predictive_recovery','seconds']:
            a=g[metric].dropna().to_numpy()
            item[metric+'_mean']=float(np.mean(a)) if len(a) else np.nan
            item[metric+'_mcse']=float(np.std(a,ddof=1)/np.sqrt(len(a))) if len(a)>1 else np.nan
        summary.append(item)
    summary=pd.DataFrame(summary);summary.to_csv(ROOT/'results/summary.csv',index=False)
    rng=np.random.default_rng(59102901);contrasts=[]
    for (s,n),g in ok.groupby(['scenario','n']):
        for comp in ['MOFA+','MOFA+ Gaussian','Rank PCA']:
            for metric in ['auc','brier','recovery']:
                w=g.pivot(index='replicate',columns='method',values=metric).dropna(subset=['Copula',comp])
                d=(w['Copula']-w[comp]).to_numpy()*( -1 if metric=='brier' else 1)
                samples=np.mean(d[rng.integers(0,len(d),(2000,len(d)))],axis=1)
                lo,hi=np.quantile(samples,[.025,.975])
                contrasts.append(dict(scenario=s,n=n,comparator=comp,metric=metric,n_pairs=len(d),difference=np.mean(d),lower=lo,upper=hi,mcse=np.std(d,ddof=1)/np.sqrt(len(d))))
    contrasts=pd.DataFrame(contrasts);contrasts.to_csv(ROOT/'results/paired_contrasts.csv',index=False)
    diagnostics=[]
    for method in ['MOFA+','MOFA+ Gaussian','Copula']:
        g=ok[ok.method==method]
        diagnostics.append(dict(method=method,successful=len(g),flagged=int((g.diagnostic_pass==False).sum()),
                                extended=int((g.extended_sampling==True).sum()) if 'extended_sampling' in g else 0,
                                median_seconds=float(g.seconds.median()),maximum_rhat=float(g.rhat_max.max()) if method=='Copula' else None))
    (ROOT/'results/diagnostics.json').write_text(json.dumps(diagnostics,indent=2))
    # Biological factor recovery is not evaluated against a fictitious
    # unchanged truth in dependence_only, whose protein pair was modified.
    fig,axes=plt.subplots(3,2,figsize=(7.4,8),sharey=True)
    for ax,s in zip(axes.flat,SCENARIOS):
        for j,m in enumerate(METHODS[:4]):
            g=summary[(summary.scenario==s)&(summary.method==m)].sort_values('n')
            ax.errorbar(np.array([0,1])+(j-1.5)*.09,g.recovery_mean,yerr=1.96*g.recovery_mcse,fmt='o-',ms=3,color=COLORS[m],label=m,capsize=2)
        ax.set(title=LABELS[s],xticks=[0,1],xticklabels=['180','360'],ylim=(0,1.02),xlabel='Training participants')
    for ax in axes[:,0]:ax.set_ylabel('Latent subspace recovery')
    axes[0,0].legend(loc='lower left',fontsize=9.5)
    fig.tight_layout();savefig(fig,'factor_recovery')
    fig,axes=plt.subplots(3,2,figsize=(7.4,8),sharey=True)
    for ax,s in zip(axes.flat,SCENARIOS):
        for j,m in enumerate(METHODS):
            g=summary[(summary.scenario==s)&(summary.method==m)].sort_values('n')
            ax.errorbar(np.array([0,1])+(j-2.5)*.065,g.auc_mean,yerr=1.96*g.auc_mcse,fmt='o-',ms=3,color=COLORS[m],label=m,capsize=2)
        ax.axhline(.5,color='grey',ls=':',lw=1)
        ax.set(title=LABELS[s],xticks=[0,1],xticklabels=['180','360'],ylim=(.4,.85),xlabel='Training participants')
    for ax in axes[:,0]:ax.set_ylabel('Independent test ROC AUC')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.07),fontsize=10)
    fig.tight_layout();savefig(fig,'prediction')
    fig,axes=plt.subplots(2,1,figsize=(7.4,6.6),sharey=False)
    for ax,n in zip(axes,[180,360]):
        for j,comp in enumerate(['MOFA+','MOFA+ Gaussian','Rank PCA']):
            g=contrasts[(contrasts.n==n)&(contrasts.metric=='auc')&(contrasts.comparator==comp)].set_index('scenario').loc[SCENARIOS]
            y=np.arange(6)+(j-1)*.18
            ax.errorbar(g.difference,y,xerr=[g.difference-g.lower,g.upper-g.difference],fmt='o',color=COLORS[comp],label=comp,capsize=2,ms=4)
        ax.axvline(0,color='grey',ls=':');ax.set(title=f'Training N = {n}',xlabel='Copula minus comparator ROC AUC',yticks=np.arange(6),yticklabels=[LABELS[s] for s in SCENARIOS])
        ax.xaxis.set_major_locator(MaxNLocator(5))
    for ax in axes:ax.invert_yaxis()
    axes[1].legend(fontsize=9,loc='best');fig.tight_layout();savefig(fig,'paired_auc')
    sample=generate(360,69102901,'missing')
    fig,axes=plt.subplots(2,2,figsize=(7.4,4.8))
    for ax,v,label in zip(axes.flat,sample['views'],['Germline dosage','Somatic indicator','RNA log2 CPM','Protein log2 abundance']):
        a=v[:,0];ax.hist(a[np.isfinite(a)],bins=15,color='#4393c3',edgecolor='white')
        ax.set(title=label,ylabel='Participants',xlabel='First feature')
    fig.tight_layout();savefig(fig,'measurements')
    # Fixed example chosen by protocol order, not the most favourable replicate.
    f=np.load(ROOT/'results/replicates/regular_180_0_fixture.npz')
    fig=plt.figure(figsize=(7.4,6));gs=fig.add_gridspec(2,3,height_ratios=[1.2,.8],width_ratios=[1,1,.055])
    axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
    for ax,key,title in zip(axes[:2],['truth','copula_correlation'],['Generating latent correlation','Estimated copula correlation']):
        im=ax.imshow(f[key],vmin=-1,vmax=1,cmap='RdBu_r');ax.set(title=title,xlabel='Molecular feature',ylabel='Molecular feature')
        for edge in [5.5,11.5,35.5]:ax.axhline(edge,color='black',lw=.5);ax.axvline(edge,color='black',lw=.5)
    fig.colorbar(im,cax=fig.add_subplot(gs[0,2]))
    tr=f['copula_trace'];axes[2].plot(tr[0,:,4],label='Chain 1',lw=.7);axes[2].plot(tr[1,:,4],label='Chain 2',lw=.7)
    axes[2].set(title='Fixed monitored correlation',xlabel='Retained draw',ylabel='Latent correlation');axes[2].legend(fontsize=9.5)
    fig.subplots_adjust(wspace=.55,hspace=.6,bottom=.12,top=.92);savefig(fig,'dependence_example')
    fig,axes=plt.subplots(1,2,figsize=(7.4,4.1))
    c=ok[ok.method=='Copula'];axes[0].hist(c.rhat_max,bins=20,color='#b2182b',alpha=.8)
    axes[0].axvline(1.1,color='black',ls=':');axes[0].set(xlabel='Maximum split Rhat\nacross 16 correlations',ylabel='Fits',title='Copula chain diagnostics')
    vals=[ok[ok.method==m].seconds for m in METHODS[:4]]
    axes[1].boxplot(vals,tick_labels=METHODS[:4],showfliers=False)
    axes[1].set(yscale='log',ylabel='Seconds per fit and prediction',title='Computation')
    for label in axes[1].get_xticklabels():label.set_rotation(45);label.set_ha('right');label.set_fontsize(9.5)
    fig.tight_layout();savefig(fig,'diagnostics')
    initial=pd.DataFrame(json.loads((ROOT/'results/dependence/records.json').read_text()))
    initial.groupby(['n','scenario','method']).reject.agg(['sum','count','mean']).to_csv(ROOT/'results/dependence_initial_summary.csv')
    dep=pd.DataFrame(json.loads((ROOT/'results/dependence_validation/records.json').read_text()));dep_summary=[]
    for (n,s,m),g in dep.groupby(['n','scenario','method']):
        k=int(g.reject.sum());lo,hi=wilson(k,len(g))
        dep_summary.append(dict(n=n,scenario=s,method=m,rejections=k,evaluations=len(g),rate=k/len(g),lower=lo,upper=hi))
    dep_summary=pd.DataFrame(dep_summary);dep_summary.to_csv(ROOT/'results/dependence_summary.csv',index=False)
    fig=plt.figure(figsize=(7.4,6));gs=fig.add_gridspec(2,2)
    axes=[fig.add_subplot(gs[0,:]),fig.add_subplot(gs[1,0]),fig.add_subplot(gs[1,1])]
    for j,m in enumerate(['Pearson log','Copula normal scores']):
        g=dep_summary[(dep_summary.scenario=='null')&(dep_summary.method==m)].sort_values('n')
        axes[0].errorbar(np.array([0,1])+(j-.5)*.08,g.rate,yerr=[g.rate-g.lower,g.upper-g.rate],fmt='o-',capsize=3,label=m)
        for ax,n in zip(axes[1:],[180,360]):
            g=dep_summary[(dep_summary.scenario!='null')&(dep_summary.method==m)&(dep_summary.n==n)].sort_values('scenario')
            ax.errorbar([.25,.4,.65],g.rate,yerr=[g.rate-g.lower,g.upper-g.rate],fmt='o-',capsize=3,label=m)
            ax.set(title=f'Detection at N = {n}',xlabel='Latent correlation in TNBC',ylabel='Rejection probability',ylim=(0,1.04))
    axes[0].set(title='Fresh null evaluation',xticks=[0,1],xticklabels=['180','360'],xlabel='Participants',ylabel='Rejection probability',ylim=(0,.12))
    axes[0].axhline(.05,color='grey',ls=':')
    axes[1].legend(fontsize=9);fig.tight_layout();savefig(fig,'dependence_test')
    # Archive computational audit with all source and output hashes.
    lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
    for name,sha in lock['sources'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha
    audit=dict(status='passed' if len(ok)==len(df) else 'failures_reported',evaluation_datasets=len(paths),method_records=len(df),failed_fits=int((df.status!='ok').sum()),
               unique_seeds=int(df.seed.nunique()),independent_test_size=200,source_lock=lock,
               cohort_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
               dependence_initial_calibration_datasets=1000,dependence_initial_evaluation_datasets=2000,
               dependence_fresh_calibration_datasets=10000,dependence_fresh_evaluation_datasets=8000)
    (ROOT/'results/audit.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps({k:v for k,v in audit.items() if k not in ['source_lock','cohort_hashes']},indent=2))


if __name__=='__main__':main()
