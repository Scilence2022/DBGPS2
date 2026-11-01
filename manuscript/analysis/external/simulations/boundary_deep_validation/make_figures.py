"""Scientific figures from frozen measurement tables; no model fitting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from deep_simulate import window_retention

R=Path(__file__).resolve().parent
BLUE='#0072B2';ORANGE='#E69F00';GREEN='#009E73';PURPLE='#CC79A7';GRAY='#777777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})


def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(R/f'{name}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.09)
    plt.close(fig)


def main():
    p=pd.read_csv(R/'policy_results.csv');b=pd.read_csv(R/'crossing_bounds_validation.csv')
    fig,gridaxes=plt.subplots(2,2,figsize=(7.2,5.8));fig.subplots_adjust(left=.10,right=.98,top=.93,bottom=.18,wspace=.33,hspace=.55)
    first=gridaxes[0,0];axs=[gridaxes[0,1],gridaxes[1,0],gridaxes[1,1]]
    metadata=[json.loads(f.read_text()) for f in (R/'runs').glob('*.input.json')]
    for geometry,color,label in [('regular',BLUE,'Regular'),('bernoulli',GREEN,'Bernoulli')]:
        cases=[m for m in metadata if m['settings']['cohort']=='poisson_theory' and m['settings']['geometry']==geometry and m['settings']['length']==25 and m['settings']['depth']==100 and m['settings']['recovery']==1]
        assert len(cases)==3
        hist=np.array([m['fragment_length_histogram'] for m in cases]).sum(axis=0);prob=hist/hist.sum()
        first.plot(np.arange(1,201),np.where(prob[1:]>0,prob[1:],np.nan),color=color,label=label)
    first.axvline(31,color=GRAY,ls=':',lw=.8);first.set_xlim(0,120);first.set_yscale('log');first.set_ylim(1e-5,1)
    first.set_xlabel('Fragment length (nt)');first.set_ylabel('Fraction of fragments');first.legend(fontsize=6.5)
    first.set_title('a  Matched cut density, different tails',loc='left',fontsize=8,fontweight='bold')
    k=21;L=45;q=.05;D=100;e=.01
    for geometry,color,label in [('regular',BLUE,'Regular'),('bernoulli',GREEN,'Bernoulli')]:
        mu=200*D/window_retention(geometry,L,q,1).sum()
        lam=mu*window_retention(geometry,L,q,k)*(1-e)**k
        axs[0].plot(np.arange(201-k)+(k-1)/2,lam,color=color,label=label)
    axs[0].axhline(2,color=GRAY,lw=.7,ls=':');axs[0].set_xlabel('Window centre (nt)');axs[0].set_ylabel('Expected origin-correct count')
    axs[0].set_title('b  Effective window support',loc='left',fontsize=8,fontweight='bold');axs[0].legend(fontsize=6.5)
    s=b[b.c==1]
    for geometry,color,label in [('regular',BLUE,'Regular'),('bernoulli',GREEN,'Bernoulli')]:
        z=s[s.geometry==geometry];mid=(z.model_lower+z.model_upper)/2
        axs[1].errorbar(100*mid,100*z.empirical,xerr=np.vstack([100*(mid-z.model_lower),100*(z.model_upper-mid)]),
                        fmt='o',ms=2.5,lw=.5,color=color,alpha=.55,label=label)
    axs[1].plot([0,100],[0,100],color=GRAY,lw=.8);axs[1].set_xlim(-2,102);axs[1].set_ylim(-2,102)
    axs[1].set_xlabel('Model probability bounds (%)');axs[1].set_ylabel('Observed complete support (%)')
    axs[1].set_title('c  Bounds using downcrossings',loc='left',fontsize=8,fontweight='bold')
    for q,color in [(1,BLUE),(.05,ORANGE)]:
        z=s[s.recovery==q];axs[2].scatter(100*z.naive_product,100*z.empirical,s=12,color=color,alpha=.6,label=f'q={q:g}')
    axs[2].plot([0,100],[0,100],color=GRAY,lw=.8);axs[2].set_xlim(-2,102);axs[2].set_ylim(-2,102)
    axs[2].set_xlabel('Product of marginals (%)');axs[2].set_ylabel('Observed complete support (%)')
    axs[2].set_title('d  Correlated windows',loc='left',fontsize=8,fontweight='bold');axs[2].legend(fontsize=6.5)
    fig.text(.10,.025,'a: L=25 nt, D=100×, q=1; vertical line k=31. b: L=45, D=100×, q=0.05, e=1%, k=21.\nc,d: c=1, three pools per condition; horizontal ranges in c are model probability bounds.',fontsize=6.6)
    save(fig,'Fig14_effective_coverage_model')

    d=p[(p.cohort=='error_structure')&(p.policy=='default')].copy()
    meta={}
    for path in (R/'runs').glob('*.input.json'):
        z=json.loads(path.read_text());meta[path.name.split('.')[0]]=z
    d['realized_error_pct']=[100*meta[case.split('_k')[0]]['realized_substitution_rate'] for case in d.case_id]
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.25));fig.subplots_adjust(left=.09,right=.98,top=.83,bottom=.25,wspace=.35)
    channels=['substitution','burst','persistent'];names=['Independent','Bursty','Persistent'];colors=[BLUE,GREEN,PURPLE]
    for j,(channel,label,color) in enumerate(zip(channels,names,colors)):
        z=d[(d.channel==channel)&(d.k==21)].groupby('depth').realized_error_pct.agg(['mean','std']).loc[[30,100]]
        axs[0].bar(np.arange(2)+(j-1)*.24,z['mean'],.22,yerr=z['std'],capsize=2,color=color,label=label)
        for k,marker,style in [(21,'o','-'),(31,'s','--')]:
            z=d[(d.channel==channel)&(d.k==k)].groupby('depth').correct_pct.agg(['mean','std']).loc[[30,100]]
            axs[1].errorbar(range(2),z['mean'],yerr=z['std'],color=color,marker=marker,mfc=color if k==21 else 'white',
                            ls=style,ms=4,capsize=2,label=label if k==21 else None)
    axs[0].set_title('a  Similar mean error fractions',loc='left',fontweight='bold');axs[0].set_ylabel('Substituted source bases (%)');axs[0].set_ylim(0,3.7)
    axs[1].set_title('b  Different recovery boundaries',loc='left',fontweight='bold');axs[1].set_ylabel('Correct strands (%)');axs[1].set_ylim(-2,105)
    for ax in axs:ax.set_xticks(range(2),['30×','100×']);ax.set_xlabel('Target retained-base depth')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,fontsize=7)
    fig.text(.09,.035,'L=80 nt, q=1, c=1, r=2.55. Mean ± s.d., three pools. b: circles k=21; open squares k=31.',fontsize=6.7)
    save(fig,'Fig15_error_structure')

    names={'default':'Fixed r=2.55','global_crc_r':'CRC-count choice','union_consensus':'Consensus union','two_support_consensus':'Two-setting support'}
    colors={'default':BLUE,'global_crc_r':ORANGE,'union_consensus':GREEN,'two_support_consensus':PURPLE}
    markers={'default':'o','global_crc_r':'s','union_consensus':'^','two_support_consensus':'D'}
    s=p[p.cohort=='scale'];fig,axs=plt.subplots(2,2,figsize=(7.2,5.8));fig.subplots_adjust(left=.10,right=.98,top=.88,bottom=.16,hspace=.5,wspace=.34)
    for row,k in enumerate([21,31]):
        z=s[s.k==k]
        for name,label in names.items():
            v=z[z.policy==name].groupby('n').agg({'correct_pct':['mean','std'],'false_accepted':['mean','std']}).loc[[1024,8192,65536]]
            for col,field in enumerate(['correct_pct','false_accepted']):
                axs[row,col].errorbar(range(3),v[(field,'mean')],yerr=v[(field,'std')] if col==0 else None,color=colors[name],
                                      marker=markers[name],mfc='none' if name in ['union_consensus','two_support_consensus'] else colors[name],
                                      ms=4,ls='--' if name=='two_support_consensus' else '-',lw=1,capsize=2,label=label)
                if col==1:
                    for x,n in enumerate([1024,8192,65536]):
                        counts=z[(z.policy==name)&(z.n==n)].sort_values('seed').false_accepted.to_numpy()
                        axs[row,col].scatter(x+np.linspace(-.05,.05,3),counts,s=9,color=colors[name],alpha=.5)
        ceiling=z[z.policy=='default'].groupby('n').coverage_ceiling_pct.mean().loc[[1024,8192,65536]]
        axs[row,0].plot(range(3),ceiling,':',color=GRAY,label='Coverage ceiling');axs[row,0].set_ylim(0,102)
        axs[row,0].set_ylabel('Correct strands (%)');axs[row,1].set_ylabel('Erroneous accepted outputs')
        from matplotlib.ticker import MaxNLocator
        axs[row,1].set_ylim(-.15,max(1.1,float(z.false_accepted.max())+.5));axs[row,1].yaxis.set_major_locator(MaxNLocator(integer=True))
        for col in [0,1]:
            axs[row,col].set_xticks(range(3),['1,024','8,192','65,536']);axs[row,col].set_xlabel('Designs per library')
            axs[row,col].set_title(f'{chr(97+row*2+col)}  k={k}: '+('recovery' if col==0 else 'errors'),loc='left',fontweight='bold')
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=3,frameon=False,fontsize=7)
    fig.text(.10,.025,'L=45 nt, D=100×, q=0.05, e=1%. Recovery: mean ± s.d. Errors: means and individual counts.\nThree independent, nested-size families.',fontsize=6.7)
    save(fig,'Fig16_blind_selection_scale')
    print('DEEP_FIGURES_COMPLETE',flush=True)


if __name__=='__main__':main()
