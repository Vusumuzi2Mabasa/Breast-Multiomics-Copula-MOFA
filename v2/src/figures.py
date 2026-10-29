"""Publication figures for the prespecified Study 3 experiments."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
ROOT=Path(__file__).resolve().parents[1]
STAT=['Partial Pearson','Partial Spearman','Normal-score difference','Product regression','Conditional copula','Oracle margins']
STAT_LABEL=['Pearson','Spearman','Normal scores','Product','Copula','Oracle']
SC=['null','abundance_only','scale_only','technical_null','missing_null','heavy_tail_null']
SC_LABEL=['No subtype change','Abundance only','Scale only','Technical variation','Missing proteins','Heavy-tail stress']
ALT=['coordination_moderate','coordination_strong','both','opposing_coordination','missing_coordination']
ALT_LABEL=['Moderate coordination','Strong coordination','Abundance and coordination','Opposing changes','Missing and coordination']
INT=['MOFA+ tuned','MOFA+ rank tuned','Gaussian Gibbs matched','Copula Gibbs matched','Rank PCA','Covariates only']
INT_LABEL=['MOFA raw','MOFA ranks','Gaussian Gibbs','Copula Gibbs','Rank PCA','Covariates']
COLORS=['#777777','#5e3c99','#2166ac','#e08214','#b2182b','#1b7837']
IC=['#2166ac','#67a9cf','#e08214','#b2182b','#5e3c99','#777777']
ISC=['abundance_only','coordination_moderate','both','missing_coordination']
IL=['Abundance only','Coordination only','Both','Missing and both']


def read(name):return pd.read_csv(ROOT/f'results/{name}.csv',keep_default_na=False,na_values=[''])


def save(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',dpi=180,bbox_inches='tight',facecolor='white')
    svg=ROOT/f'figures/{name}.svg'
    fig.savefig(svg,bbox_inches='tight',facecolor='white');plt.close(fig)
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')


def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    d=read('dependence_summary');i=read('integration_summary');contr=read('dependence_contrasts');raw=read('integration_evaluation');tune=read('integration_tuning')
    fig,ax=plt.subplots(figsize=(7.4,4.9));ax.set_xlim(0,10);ax.set_ylim(0,7);ax.axis('off')
    for x,y,text in [(1.4,5.8,'Germline dosages\n8 synthetic variants'),(7.9,5.8,'Somatic alterations\n8 binary indicators'),(2.1,3.2,'RNA\n16 matched genes'),(7.5,3.2,'Protein\n16 matched proteins')]:
        ax.text(x,y,text,ha='center',va='center',fontsize=12)
    for start,end in [((1.7,5.2),(2.1,3.8)),((7.3,5.2),(2.7,3.8)),((2,5.2),(7,3.8)),((7.9,5.2),(7.5,3.8))]:ax.annotate('',xy=end,xytext=start,arrowprops=dict(arrowstyle='->',color='#555555',lw=1.4))
    ax.annotate('',xy=(6.3,3.2),xytext=(3.4,3.2),arrowprops=dict(arrowstyle='<->',lw=2,color='#b2182b'))
    ax.text(4.8,3.8,'Conditional coordination',ha='center',color='#b2182b')
    ax.text(4.8,2.4,'Four prespecified modules of four gene–protein pairs',ha='center')
    ax.text(4.8,1.5,'Module 1 changes dependence in alternatives\nModules 2–4 remain unchanged',ha='center',linespacing=1.3)
    ax.text(4.8,.5,'Age, cell composition, treatment, country and batch\naffect measurements and/or subtype sampling',ha='center',fontsize=10)
    save(fig,'biological_design')
    fig,axes=plt.subplots(2,2,figsize=(7.4,6.8),sharex=True)
    for row,n in enumerate([240,480]):
        for col,metric in enumerate(['false_module','false_gene']):
            ax=axes[row,col]
            for j,(m,label,color) in enumerate(zip(STAT,STAT_LABEL,COLORS)):
                g=d[(d.n==n)&(d.method==m)].set_index('scenario').loc[SC]
                ax.errorbar(np.arange(6)+(j-2.5)*.055,g[metric],yerr=[g[metric]-g[metric+'_lower'],g[metric+'_upper']-g[metric]],fmt='o',ms=3,capsize=1.5,color=color,label=label)
            ax.axhline(.05,color='black',ls=':',lw=1);ax.set(title=f'{"Module" if col==0 else "Gene"} family, N={n}',ylabel='Any false discovery',ylim=(0,max(.15,float(d[(d.n==n)&d.scenario.isin(SC)][metric+'_upper'].max())*1.08)))
            ax.set_xticks(np.arange(6));ax.set_xticklabels(['Null','Mean','Scale','Tech','Missing','t stress'],rotation=40,ha='right')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.07),fontsize=10)
    fig.tight_layout();save(fig,'null_control')
    fig,axes=plt.subplots(2,1,figsize=(7.4,6.4))
    for ax,n in zip(axes,[240,480]):
        for j,(m,label,color) in enumerate(zip(STAT,STAT_LABEL,COLORS)):
            g=d[(d.n==n)&(d.method==m)].set_index('scenario').loc[ALT]
            ax.errorbar(g.module_detection,np.arange(5)+(j-2.5)*.10,xerr=[g.module_detection-g.module_detection_lower,g.module_detection_upper-g.module_detection],fmt='o',ms=3,capsize=1.5,color=color,label=label)
        ax.set(title=f'Altered module detection, N={n}',yticks=np.arange(5),yticklabels=ALT_LABEL,xlabel='Detection probability',xlim=(-.03,1.03));ax.invert_yaxis()
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.07),fontsize=10)
    fig.tight_layout();save(fig,'module_detection')
    fig,axes=plt.subplots(2,1,figsize=(7.4,6.1))
    for ax,n in zip(axes,[240,480]):
        for j,(comp,color,label) in enumerate([('Product regression','#e08214','Product'),('Normal-score difference','#2166ac','Normal scores'),('Partial Spearman','#5e3c99','Spearman')]):
            g=contr[(contr.n==n)&(contr.comparator==comp)].set_index('scenario').loc[ALT]
            ax.errorbar(g.difference,np.arange(5)+(j-1)*.15,xerr=[np.maximum(g.difference-g.lower,0),np.maximum(g.upper-g.difference,0)],fmt='o',capsize=2,ms=3,color=color,label=label)
        ax.axvline(0,color='grey',ls=':');ax.set(title=f'Paired detection difference, N={n}',yticks=np.arange(5),yticklabels=ALT_LABEL,xlabel='Copula minus comparator');ax.invert_yaxis();ax.xaxis.set_major_locator(MaxNLocator(5))
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.04),fontsize=9)
    fig.tight_layout();save(fig,'paired_detection')
    fig,axes=plt.subplots(2,2,figsize=(7.4,5.7))
    for ax,method in zip(axes[0],['MOFA+ tuned','MOFA+ rank tuned']):
        g=raw[raw.method==method];ax.hist(g.k,bins=[1,3,5,8,12],rwidth=.7,color='#4393c3');ax.set(title=method,xlabel='Selected factor budget',ylabel='Replicates',xticks=[2,4,6,10])
    for method in ['Gaussian Gibbs matched','Copula Gibbs matched']:
        g=raw[raw.method==method];axes[1,0].hist(g.rhat_max,bins=np.linspace(1, max(1.12,raw.rhat_max.max()+.01),18),alpha=.55,label=method.split()[0])
    axes[1,0].axvline(1.1,color='black',ls=':');axes[1,0].set(xlabel='Maximum monitored split Rhat',ylabel='Replicates',title='Matched sampler diagnostics');axes[1,0].legend(fontsize=9)
    methods=INT[:5];axes[1,1].boxplot([raw[raw.method==m].seconds for m in methods],tick_labels=['MOFA','Rank\nMOFA','Gauss\nGibbs','Copula\nGibbs','PCA'],showfliers=False)
    axes[1,1].set(yscale='log',ylabel='Seconds per workflow',title='Selection and fitting cost');axes[1,1].tick_params(axis='x',labelsize=9)
    fig.tight_layout();save(fig,'integration_diagnostics')
    fig,axes=plt.subplots(2,2,figsize=(7.4,6.5))
    for ax,metric,title in zip(axes.flat,['auc','capture','protein_mse','balanced_accuracy'],['Subtype discrimination','Generating module capture','Missing protein-view prediction','Supervised comparison']):
        if metric=='balanced_accuracy':
            di=read('diablo_evaluation');dg=di[di.n==480].groupby('scenario').balanced_accuracy.agg(['mean','std','count']).loc[ISC]
            ax.errorbar(np.arange(4)-.15,dg['mean'],yerr=1.96*dg['std']/np.sqrt(dg['count']),fmt='s',capsize=2,color='black',label='DIABLO')
            plot_methods=['MOFA+ tuned','Copula Gibbs matched','Covariates only']
        else:plot_methods=INT if metric!='capture' else INT[:5]
        for j,m in enumerate(plot_methods):
            g=i[(i.n==480)&(i.method==m)].set_index('scenario').loc[ISC]
            ax.errorbar(np.arange(4)+(j-(len(plot_methods)-1)/2)*.045,g[metric],yerr=1.96*g[metric+'_mcse'],fmt='o',ms=3,capsize=1.5,color=IC[INT.index(m)],label=INT_LABEL[INT.index(m)])
        ax.set(title=title,xticks=np.arange(4),xticklabels=['Mean','Coord','Both','Missing'],ylabel={'auc':'ROC AUC','capture':'Subspace capture','protein_mse':'Standardised mean squared error','balanced_accuracy':'Balanced accuracy'}[metric])
        if metric in ['auc','balanced_accuracy']:ax.axhline(.5,color='grey',ls=':')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.07),fontsize=9.5)
    axes[1,1].legend(fontsize=8,loc='best');fig.tight_layout();save(fig,'integration_results')
    transport=read('transport_summary')
    fig,axes=plt.subplots(1,2,figsize=(7.4,4.5),sharey=True)
    for ax,n in zip(axes,[240,480]):
        for j,m in enumerate(INT[:5]):
            g=transport[(transport.n==n)&(transport.method==m)].set_index('scenario').loc[ISC]
            ax.errorbar(np.arange(4)+(j-2)*.04,g.difference,yerr=[g.difference-g.lower,g.upper-g.difference],fmt='o',capsize=1.5,ms=3,color=IC[j],label=INT_LABEL[j])
        ax.axhline(0,color='grey',ls=':');ax.set(title=f'Training N={n}',xticks=np.arange(4),xticklabels=['Mean','Coord','Both','Missing'],ylabel='Unseen minus seen-country AUC')
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',ncol=3,bbox_to_anchor=(.5,1.1),fontsize=9.5)
    fig.tight_layout();save(fig,'country_transport')
    print('Seven publication figures generated.')


if __name__=='__main__':main()
