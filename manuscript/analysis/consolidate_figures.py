"""Seven main figures from preserved measurements; no new reconstruction runs.

Figures 2-6 and 8. Reads the snapshots in ../data/source_data and writes to ../figures;
Figure 1 comes from fig1_concept.py and Figure 7 from fig7_ht176.py.
"""
from pathlib import Path
import gzip,hashlib,json,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator,FuncFormatter

M=Path(__file__).resolve().parents[1];R=M/'figures';O=M
DATA=M/'analysis/data/source_data'
T=Path('/Users/song/dbgps2_ht/results')
TD=M/'analysis/data'
B=O/'exploration/2026-10-20/boundary_deep_validation'
EXT=O/'figures/2026-10-19/DBGPS2_external_figures/data'
BENCH=O/'benchmarks/2026-10-19'
MAN=O/'docx/DBGPS2_manuscript_recovery_speed_restructured.docx'
manifest={}
def snapshot(p,label=None):
    p=Path(p);label=label or p.name;d=DATA/label
    if not d.exists():shutil.copyfile(p,d)
    manifest[label]={'source':str(p),'bytes':d.stat().st_size,'sha256':hashlib.sha256(d.read_bytes()).hexdigest()}
    return d
def csv(p,label=None,**kwargs):return pd.read_csv(snapshot(p,label),**kwargs)
def js(p,label=None):return json.loads(snapshot(p,label).read_text())

BLUE='#0072B2';ORANGE='#D46A36';GREEN='#009E73';PURPLE='#CC79A7';GRAY='#7B8389';LIGHT='#D9DEE2';RED='#B44D43';INK='#25343E'
SAMPLES=['HT_con','HT1','HT2','HT3','HT4','HT5'];LABELS=['Control','HT1','HT2','HT3','HT4','HT5']
HC=dict(zip(SAMPLES,['#3b3b3b','#f2a07b','#e57850','#d1502f','#a8321c','#771b10']))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.2,'axes.titlesize':8,'axes.labelsize':7.2,
    'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.2,'axes.spines.top':False,
    'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.1,'pdf.fonttype':42,'svg.fonttype':'none'})
def title(ax,letter,label):ax.set_title(f'{letter}  {label}',loc='left',fontweight='bold',pad=9)
def save(fig,n,name):
    for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig{n}_{name}.{ext}',dpi=400,bbox_inches='tight',pad_inches=.10)
    plt.close(fig)
def grid(rows=3,cols=2,heights=None):
    fig,axs=plt.subplots(rows,cols,figsize=(7.2,6.25 if rows==3 else 5.1),gridspec_kw={'height_ratios':heights} if heights else {})
    fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.10,hspace=.85,wspace=.43)
    return fig,np.asarray(axs).ravel()

# Figure 6 combines the three simulation figures into six mechanistic panels.
p=csv(B/'policy_results.csv');bound=csv(B/'crossing_bounds_validation.csv')
fig,ax=grid()
for geom,col,lab in [('regular',BLUE,'Regular'),('bernoulli',GREEN,'Bernoulli')]:
    selected=[]
    for f in (B/'runs').glob('*.input.json'):
        m=json.loads(f.read_text());s=m['settings']
        if s['cohort']=='poisson_theory' and s['geometry']==geom and s['length']==25 and s['depth']==100 and s['recovery']==1:
            selected.append(js(f,'simulation_'+f.name))
    if not selected: # Portable redraw from the six archived metadata snapshots.
        for f in DATA.glob('simulation_*.input.json'):
            m=json.loads(f.read_text())
            if m['settings']['geometry']==geom:selected.append(m)
    assert len(selected)==3
    h=np.array([m['fragment_length_histogram'] for m in selected]).sum(0);prob=h/h.sum()
    ax[0].plot(np.arange(1,201),np.where(prob[1:]>0,prob[1:],np.nan),color=col,label=lab)
ax[0].axvline(31,color=GRAY,ls=':',lw=.8);ax[0].set(xlim=(0,120),ylim=(1e-5,1),yscale='log',xlabel='Fragment length (nt)',ylabel='Fraction of fragments');ax[0].legend();title(ax[0],'a','Fragment tails preserve long windows')
s=bound[bound.c==1]
for geom,col in [('regular',BLUE),('bernoulli',GREEN)]:
    z=s[s.geometry==geom];mid=(z.model_lower+z.model_upper)/2
    ax[1].errorbar(100*mid,100*z.empirical,xerr=np.vstack([100*(mid-z.model_lower),100*(z.model_upper-mid)]),fmt='o',ms=2.4,lw=.4,color=col,alpha=.55)
ax[1].plot([0,100],[0,100],color=GRAY,lw=.7);ax[1].set(xlim=(-2,102),ylim=(-2,102),xlabel='Dependent window probability bounds (%)',ylabel='Observed complete support (%)');title(ax[1],'b','Bounds account for shared support')
for q,col in [(1,BLUE),(.05,'#E69F00')]:
    z=s[s.recovery==q];ax[2].scatter(100*z.naive_product,100*z.empirical,s=10,color=col,alpha=.65,label=f'q = {q:g}')
ax[2].plot([0,100],[0,100],color=GRAY,lw=.7);ax[2].set(xlim=(-2,102),ylim=(-2,102),xlabel='Product of marginal probabilities (%)',ylabel='Observed complete support (%)');ax[2].legend();title(ax[2],'c','Independent windows misstate availability')
err=p[(p.cohort=='error_structure')&(p.policy=='default')]
for channel,col,lab in [('substitution',BLUE,'Independent'),('burst',GREEN,'Bursty'),('persistent',PURPLE,'Persistent')]:
    for k,marker,ls in [(21,'o','-'),(31,'s','--')]:
        z=err[(err.channel==channel)&(err.k==k)].groupby('depth').correct_pct.agg(['mean','std']).loc[[30,100]]
        ax[3].errorbar(range(2),z['mean'],yerr=z['std'],color=col,marker=marker,mfc=col if k==21 else 'white',ls=ls,ms=3.5,capsize=2,label=lab if k==21 else None)
ax[3].set(xticks=[0,1],xticklabels=['30×','100×'],ylim=(-2,105),xlabel='Target retained-base depth',ylabel='Correct strands (%)');ax[3].legend(loc='center right',fontsize=6);title(ax[3],'d','Error correlation shifts recovery')
policies={'default':('Fixed r = 2.55',BLUE,'o'), 'global_crc_r':('CRC-count choice','#E69F00','s'),'union_consensus':('Consensus union',GREEN,'^'),'two_support_consensus':('Two-setting support',PURPLE,'D')}
sc=p[(p.cohort=='scale')&(p.k==21)]
for name,(lab,col,mk) in policies.items():
    z=sc[sc.policy==name];v=z.groupby('n').agg({'correct_pct':['mean','std'],'false_accepted':['mean']}).loc[[1024,8192,65536]]
    ax[4].errorbar(range(3),v[('correct_pct','mean')],yerr=v[('correct_pct','std')],color=col,marker=mk,ms=3.5,capsize=2,label=lab)
    ax[5].plot(range(3),v[('false_accepted','mean')],color=col,marker=mk,ms=3.5)
    for x,n in enumerate([1024,8192,65536]):
        vals=z[z.n==n].sort_values('seed').false_accepted.to_numpy();ax[5].scatter(x+np.linspace(-.05,.05,3),vals,s=8,color=col,alpha=.45)
ceiling=sc[sc.policy=='default'].groupby('n').coverage_ceiling_pct.mean().loc[[1024,8192,65536]]
ax[4].plot(range(3),ceiling,':',color=GRAY,label='Output-support ceiling');ax[4].set(ylim=(0,103),ylabel='Correct strands (%)');ax[4].legend(loc='lower left',fontsize=5.9,ncol=2)
ax[5].set(ylim=(-.15,5.6),ylabel='Erroneous accepted outputs');ax[5].yaxis.set_major_locator(MaxNLocator(integer=True))
for a in ax[4:]:a.set(xticks=range(3),xticklabels=['1,024','8,192','65,536'],xlabel='Designs per library')
title(ax[4],'e','Blind tuning increases correct yield');title(ax[5],'f','CRC agreement retains false paths')
save(fig,6,'simulation_boundaries')

# Figure 2 connects physical fragmentation, graph interference and rescued strands.
fig,ax=grid()
for samp,lab in zip(SAMPLES,LABELS):
    fr=csv(T/f'qc/{samp}.frag.tsv',f'thermal_{samp}.frag.tsv',sep='\t').to_numpy();h=fr[:,1]
    ax[0].step(fr[:,0],np.cumsum(h)/h.sum()*100,where='post',color=HC[samp],label=lab)
    pos=sum(csv(T/f'qc/{samp}_{mate}.pos.tsv',f'thermal_{samp}_{mate}.pos.tsv',sep='\t').to_numpy(dtype=float) for mate in [1,2]);pos[:,0]/=2
    ax[1].plot(pos[:,0],pos[:,1]/pos[18:182,1].max(),color=HC[samp])
ax[0].set(xlim=(0,205),ylim=(0,102),xlabel='Aligned fragment length (nt)',ylabel='Aligned pairs ≤ length (%)');ax[0].legend(ncol=2,loc='upper left',fontsize=6);title(ax[0],'a','Aging leaves short recovered fragments')
ax[1].axvspan(0,18,color=LIGHT,alpha=.6,lw=0);ax[1].axvspan(182,200,color=LIGHT,alpha=.6,lw=0);ax[1].set(xlim=(0,200),ylim=(0,1.07),xlabel='Position along 200-nt strand',ylabel='Relative read depth');title(ax[1],'b','Central coverage is depleted')
for samp in ['HT_con','HT5']:
    spec=csv(T/f'full/{samp}.spectrum.tsv',f'thermal_{samp}.spectrum.tsv',sep='\t',comment='#')
    for field,ls,label in [('ref','-', 'true'),('noise','--','noise')]:ax[2].plot(spec['cov'],spec[field],color=HC[samp],ls=ls,label=f'{dict(zip(SAMPLES,LABELS))[samp]} {label}')
    rh=csv(T/f'full/{samp}.ratiohist.tsv',f'thermal_{samp}.ratiohist.tsv',sep='\t')
    for field,ls in [('true','-'),('alt_ge2','--')]:ax[3].plot(rh.log2_ratio_lo+.025,rh[field]/rh[field].sum()*100,color=HC[samp],ls=ls)
ax[2].set(xscale='log',yscale='log',xlim=(1,3000),ylim=(1,2e9),xlabel='k-mer coverage',ylabel='Distinct k-mers');ax[2].legend(fontsize=5.9);title(ax[2],'c','Absolute counts overlap with noise')
for v in [-np.log2(2.55),np.log2(2.55)]:ax[3].axvline(v,color=BLUE,lw=.7)
ax[3].set(yscale='log',xlim=(-10,3),ylim=(1e-4,3000),xlabel='log₂(child count / parent count)',ylabel='Transitions (%)');ax[3].text(.97,.97,'Solid: true transitions\nDashed: eligible noise',transform=ax[3].transAxes,ha='right',va='top',fontsize=6);title(ax[3],'d','Local ratios separate true transitions')
stats=csv(T/'full/HT5.strandstats.tsv','thermal_HT5.strandstats.tsv.gz',sep='\t')
trace=csv(T/'full/HT5.trace.tsv','thermal_HT5.trace.tsv',sep='\t');ss=stats.set_index('id');cand=[i for i in trace.id.unique() if ss.loc[i,'min_cov']>=3 and ss.loc[i,'n_alt_ge2']>=3];ex=cand[0] if cand else trace.id.iloc[0];tr=trace[trace.id==ex]
ax[4].plot(tr.step,tr.true_cov,color=INK,label='True path');ax[4].fill_between(tr.step,tr.parent_cov/2.55,tr.parent_cov*2.55,color=BLUE,alpha=.15,label='r = 2.55 window')
for i,c in enumerate(['alt1','alt2','alt3']):
    m=tr[c]>0;ax[4].scatter(tr.step[m],tr[c][m],s=5,color=ORANGE,label='Noise' if i==0 else None)
ax[4].set(yscale='log',ylim=(.7,3000),xlabel='Extension step',ylabel='k-mer coverage');ax[4].legend(fontsize=5.8);title(ax[4],'e','Local windows follow a weak true path')
out=snapshot(T/'full/HT5.outcomes.tsv','thermal_HT5.outcomes.tsv.gz');hdr=gzip.open(out,'rt').readline().strip().split('\t')[1:];oc=pd.read_csv(out,sep='\t',skiprows=1,header=None,names=['id','o'])
g=csv(T/'full/HT5.grid.tsv','thermal_HT5.grid.tsv',sep='\t');d=g[g.r==0].sort_values('correct',ascending=False).iloc[0]
o1=oc.o.str[hdr.index(f'c{int(d.c)}_r0')].to_numpy();o2=oc.o.str[hdr.index('c1_r2.55')].to_numpy();mc=ss.loc[oc.id,'min_cov'].to_numpy()
for lab,mask,col in [('Both',(o1=='S')&(o2=='S'),'#7BB6DF'),('DBGPS2 only',(o1!='S')&(o2=='S'),BLUE),('Neither',(o1!='S')&(o2!='S'),GRAY)]:
    vals=np.sort(np.maximum(mc[mask],.5));ax[5].plot(vals,np.arange(1,len(vals)+1)/len(vals)*100,color=col,label=lab)
ax[5].axvline(float(d.c)+.5,color=ORANGE,ls='--',lw=.8);ax[5].set(xscale='log',xlim=(.5,300),ylim=(0,102),xlabel='Minimum reference-path count',ylabel='Strands in group (cumulative %)');ax[5].legend(loc='lower right',fontsize=5.8);title(ax[5],'f','DBGPS2 rescues weakly covered strands')
save(fig,2,'degradation_mechanism')

# Figure 3 combines recovery, parameter choice and sequencing requirement.
grids={s:csv(T/f'full/{s}.grid.tsv',f'thermal_{s}.grid.tsv',sep='\t') for s in SAMPLES}
kt=csv(TD/'ksize_sweep.csv','thermal_ksize_sweep.csv');dep=csv(TD/'depth_series.csv','thermal_depth_series.csv');need=csv(TD/'depth_needed_95.csv','thermal_depth_needed_95.csv')
fig,ax=grid();x=np.arange(6);w=.36
for j,(ratio,col,label) in enumerate([(0,ORANGE,'DBGPS rule'),(2.55,BLUE,'DBGPS2')]):
    y=[grids[s][grids[s].r==ratio].correct.max()/2100 for s in SAMPLES];ax[0].bar(x+(j-.5)*w,y,w,color=col,label=label)
ax[0].set(xticks=x,xticklabels=LABELS,ylim=(0,125),yticks=[0,20,40,60,80,100],ylabel='Exact strands (%)');ax[0].legend(loc='upper right',ncol=2,fontsize=6,frameon=False);title(ax[0],'a','Recovery advantage grows with damage')
g=grids['HT5'];rs=sorted(g.r.unique(),key=lambda v:(v==0,v));cs=sorted(g.c.unique());mat=np.array([[g[(g.r==r)&(g.c==c)].correct.iloc[0]/2100 for c in cs] for r in rs])
im=ax[1].imshow(mat,aspect='auto',vmin=0,vmax=100,cmap='Blues',origin='lower');ax[1].set(xticks=range(len(cs)),xticklabels=cs,yticks=range(len(rs)),yticklabels=['Off' if r==0 else f'{r:g}' for r in rs],xlabel='Absolute cutoff c',ylabel='Ratio r');ax[1].tick_params(axis='x',labelrotation=90,labelsize=5.5);fig.colorbar(im,ax=ax[1],fraction=.045,pad=.025,label='Exact (%)');title(ax[1],'b','HT5 has a broad useful ratio range')
cats=[('correct','Exact',BLUE),('no_path','Broken path',LIGHT),('overflow','Overflow',ORANGE),('no_valid','No valid path','#E69F00'),('ambiguous','Ambiguous','#8272B2'),('wrong','Wrong',RED)]
for i,s in enumerate(SAMPLES):
    for j,r in enumerate([0,2.55]):
        z=grids[s][grids[s].r==r];row=z.loc[z.correct.idxmax()];base=0
        for field,label,col in cats:
            v=row[field]/2100;ax[2].bar(i*2.5+j,v,.85,bottom=base,color=col,label=label if i==0 and j==0 else None);base+=v
ax[2].set(xticks=x*2.5+.5,xticklabels=LABELS,ylim=(0,102),ylabel='Designed strands (%)');ax[2].legend(loc='lower center',bbox_to_anchor=(.5,1.20),ncol=3,fontsize=5.6,handlelength=.8,frameon=False);title(ax[2],'c','Local pruning suppresses path overflow')
for s in SAMPLES:
    z=kt[kt['sample']==s].sort_values('k');ax[3].plot(z.k,z.dbgps2/2100,color=HC[s],marker='o',ms=2.5,label=dict(zip(SAMPLES,LABELS))[s])
ax[3].set(xticks=[17,19,21,23,25,27,31],xlabel='k-mer size',ylabel='Exact strands by DBGPS2 (%)');ax[3].legend(ncol=3,fontsize=5.6,loc='lower center');title(ax[3],'d','Shorter windows improve continuity')
for r,k,col,ls,label in [(0,31,ORANGE,'-','DBGPS k31'),(2.55,31,BLUE,'-','DBGPS2 k31'),(2.55,21,BLUE,'--','DBGPS2 k21')]:
    z=dep[(dep['sample']=='HT5')&(dep.r==r)&(dep.k==k)].sort_values('pct');reads=z.pct*2*58.64*1e6/100/210000
    ax[4].plot(reads,z.best/2100,color=col,ls=ls,marker='o',ms=2.2,label=label)
ax[4].set(xscale='log',ylim=(0,102),xlabel='Reads per designed strand',ylabel='Exact strands (%)');ax[4].legend(loc='upper left',fontsize=5.7);title(ax[4],'e','HT5 recovery increases with read depth')
for j,(lab,col) in enumerate([('DBGPS k31',ORANGE),('DBGPS2 k31',BLUE),('DBGPS2 k21','#7BB6DF')]):
    vals=np.array([need[(need['sample']==name)&(need.method==lab)].reads_per_strand_for_95pct.iloc[0] for name in LABELS]);ax[5].bar(x+(j-1)*.25,np.nan_to_num(vals),.24,color=col,label=lab)
    for xx,v in zip(x+(j-1)*.25,vals):
        if np.isnan(v):ax[5].text(xx,1.2,'n.r.',rotation=90,ha='center',fontsize=5)
ax[5].set(yscale='log',ylim=(1,8000),xticks=x,xticklabels=LABELS,ylabel='Reads per strand for 95% recovery');ax[5].legend(ncol=3,fontsize=5.6,loc='upper center',columnspacing=.8,handlelength=1,frameon=False);title(ax[5],'f','Lower sequencing requirement')
save(fig,3,'degraded_recovery')

# Figure 4 separates strand yield, the per-block outer-code requirement and the weak block of File 3.
d=js(EXT/'figure_data.json');photo_bench=csv(BENCH/'reconstruction_speed/Table5_reconstruction_speed.csv')
blk=csv(DATA/'photolith_file3_blocks.csv');pos=csv(DATA/'photolith_file3_position.csv');NEED=10977;PALE='#7BB6DF'
fig=plt.figure(figsize=(7.2,6.25));gs=fig.add_gridspec(3,6,left=.085,right=.985,top=.95,bottom=.085,hspace=.85,wspace=1.75)
ax=[fig.add_subplot(g) for g in (gs[0,:3],gs[0,3:],gs[1,:3],gs[1,3:],gs[2,:2],gs[2,2:4],gs[2,4:])]
correct=np.array([z['correct']/z['designed']*100 for z in d['photolithographic']]);wrong=np.array([z['wrong']/z['designed']*100 for z in d['photolithographic']]);xx=np.arange(3)
ax[0].bar(xx,correct,.58,color=BLUE,label='Exact');ax[0].bar(xx,wrong,.58,bottom=correct,color=RED,label='Wrong');ax[0].bar(xx,100-correct-wrong,.58,bottom=correct+wrong,color=LIGHT,label='Missing')
for i,z in enumerate(d['photolithographic']):ax[0].text(i,correct[i]/2,f'{correct[i]:.2f}%',ha='center',color='white',fontsize=7);ax[0].text(i,103,f"{z['wrong']} wrong",ha='center',fontsize=6)
ax[0].set(xticks=xx,xticklabels=['File 1\nfile recovered','File 2\nfile recovered','File 3\nfile incomplete'],ylim=(0,130),yticks=[0,25,50,75,100],ylabel='Designed strands (%)');ax[0].legend(ncol=3,loc='upper center',bbox_to_anchor=(.5,1.04),fontsize=5.8,frameon=False);title(ax[0],'a','Pooled reconstruction at high error rates')
for j,(method,col,label) in enumerate([('cluster_msa_vote',GRAY,'MSA'),('dbgps_rule',ORANGE,'DBGPS rule'),('dbgps2',BLUE,'DBGPS2')]):
    z=photo_bench[photo_bench.method==method].sort_values('dataset');ax[1].bar(np.arange(2)+(j-1)*.24,z.exact_recovery*100,.23,color=col,label=label)
ax[1].set(xticks=[0,1],xticklabels=['File 1','File 2'],ylim=(0,125),yticks=[0,25,50,75,100],ylabel='Exact strands (%)');ax[1].legend(ncol=3,loc='upper center',fontsize=5.8,frameon=False);title(ax[1],'b','Matched graph and MSA reconstruction')
total=[12*NEED,int(blk.pooled_correct.sum()),int(blk.routed_correct.sum())];ys=[2,1,0]
ax[2].barh(ys,total,.6,color=[GRAY,PALE,BLUE]);ax[2].axvline(12*NEED,color=RED,ls='--',lw=.8)
for y,v in zip(ys,total):ax[2].text(v+2500,y,f'{v:,}',va='center',fontsize=6.5)
ax[2].set(yticks=ys,yticklabels=['Code requires','Pooled output','Routed output'],xlim=(0,215000),xlabel='File 3 exact strands (all blocks)');ax[2].xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));title(ax[2],'c','Total yield exceeds the requirement')
bx=blk.block.to_numpy();ax[3].bar(bx-.2,blk.pooled_correct,.38,color=PALE,label='Pooled');ax[3].bar(bx+.2,blk.routed_correct,.38,color=[BLUE]*11+[RED],label='Routed');ax[3].axhline(NEED,color=RED,ls='--',lw=.8,label='Required per block (10,977)')
ax[3].text(12.2,NEED+450,'−1,004',ha='center',fontsize=5.8,color=RED)
ax[3].set(xticks=bx,xlim=(.2,13),ylim=(0,21000),xlabel='File 3 outer-code block',ylabel='Exact strands');ax[3].yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));ax[3].legend(ncol=3,loc='upper center',fontsize=5.6,frameon=False,columnspacing=.9,handlelength=1.2);title(ax[3],'d','Redundancy cannot cross blocks')
ax[4].bar(bx,blk.assigned_reads/1e6,.7,color=[GRAY]*11+[RED]);ax[4].set(xticks=[1,4,8,12],ylim=(0,5.4),xlabel='File 3 block',ylabel='Assigned reads (million)');title(ax[4],'e','Few reads in block 12')
ax[4].text(12,float(blk.assigned_reads.iloc[11])/1e6+.15,'1.68',ha='center',fontsize=6,color=RED)
px=(pos.block-1+(pos.tenth-.5)/10).to_numpy();py=pos.percent.to_numpy();ax[5].axvspan(11,12,color=RED,alpha=.1,lw=0)
for b in range(1,12):ax[5].axvline(b,color=LIGHT,lw=.4)
ax[5].plot(px[:111],py[:111],color=BLUE,lw=.9);ax[5].plot(px[110:],py[110:],color=RED,lw=1.2)
ax[5].set(xticks=[.5,3.5,7.5,11.5],xticklabels=['1','4','8','12'],xlim=(0,12),ylim=(0,100),xlabel='File 3 block, in index order',ylabel='Correct candidates (%)');title(ax[5],'f','Decline along the index')
show=blk[blk.block.isin([1,6,11,12])];base=np.zeros(4);lab=[str(b) for b in show.block]
for field,col,name in [('absent_kmer',GRAY,'True k-mer absent'),('singleton_only',LIGHT,'Singleton support'),('overflow',ORANGE,'Path overflow'),('ratio_rejected',BLUE,'Ratio rejected')]:
    ax[6].bar(lab,show[field],.6,bottom=base,color=col,label=name);base+=show[field].to_numpy()
ax[6].set(ylim=(0,17500),yticks=[0,5000,10000],xlabel='Pooled blocks below bound',ylabel='Indices lacking a candidate');ax[6].yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));ax[6].legend(ncol=2,loc='upper center',fontsize=4.8,frameon=False,columnspacing=.6,handlelength=.9,handletextpad=.3);title(ax[6],'g','Block 12 lacks evidence')
save(fig,4,'high_error_reconstruction')

# Figure 5 reports exact DNA Fountain and Organick strands; file recovery is a separate endpoint.
fig=plt.figure(figsize=(7.2,4.5));gs=fig.add_gridspec(2,2,height_ratios=[1,.45],hspace=.65,wspace=.45);fig.subplots_adjust(left=.10,right=.98,top=.93,bottom=.06)
a=fig.add_subplot(gs[0,0]);bax=fig.add_subplot(gs[0,1]);cax=fig.add_subplot(gs[1,:]);xx=np.arange(2)
for shift,key,col,label in [(-.17,'merge',ORANGE,'Read merging'),(.17,'DBGPS2',BLUE,'DBGPS2')]:
    vals=[d['fountain'][s][key] for s in ['full','one_tile']];a.bar(xx+shift,vals,.32,color=col,label=label)
    for x,y in zip(xx+shift,vals):a.text(x,y+900,f'{y:,}',ha='center',fontsize=6.5)
a.set(xticks=xx,xticklabels=['Full depth','One tile'],ylim=(0,90000),ylabel='Distinct strands matching a design');a.yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1000:.0f}k'));a.legend(ncol=2,loc='upper center',fontsize=6);title(a,'a','DNA Fountain exact strands')
for row,col,offset in zip(d['organick'],[LIGHT,INK,ORANGE,BLUE,GREEN,'#8272B2'],[(-10,6),(6,-13),(-13,-2),(-12,-1),(6,-2),(6,5)]):
    bax.scatter(row['recovery_pct'],row['wrong'],color=col,s=30,edgecolors='white',lw=.5);bax.annotate(row['id'],(row['recovery_pct'],row['wrong']),xytext=offset,textcoords='offset points',fontsize=7,fontweight='bold')
bax.set(xlim=(95,98.5),ylim=(-8,190),xlabel='Exact 110-nt bodies / all designs (%)',ylabel='Incorrect outputs');title(bax,'b','Organick recovery and precision')
cax.axis('off');title(cax,'c','File validation remains a separate endpoint')
cax.text(0,.65,'DNA Fountain',fontweight='bold',fontsize=8);cax.text(.30,.65,'71,990 exact strands; 2,146,816-byte archive, published MD5 matched',color=GREEN,fontsize=7.5)
cax.text(0,.18,'Organick File 16',fontweight='bold',fontsize=8);cax.text(.30,.18,'Outer-code decoding and original-file recovery not tested',color=INK,fontsize=7.5)
save(fig,5,'format_transfer')

# Figure 8 replaces six partly overlapping performance figures with four panels.
rc=csv(BENCH/'resource_control/all_configurations_summary.csv');mem=csv(BENCH/'speed_extension/Table6_configuration_comparison.csv')
ds=['fountain','nbt_file16','song_ht176'];dl=['DNA\nFountain','Organick\nFile 16','HT176']
cases=['dbgps_rule','dbgps2','dbgps2_alt_c','dbgps2_B8','dbgps2_B16','dbgps2_F32','dbgps2_F34','dbgps2_B12_F34']
names=['DBGPS rule','DBGPS2 primary c','DBGPS2 alternate c','Query B8','Query B16','Counting F32','Counting F34','B12 and F34']
get=lambda frame,w,ds,c:frame[(frame.workers==w)&(frame.dataset==ds)&(frame.case==c)].iloc[0]
fig,ax=grid(2,2)
for a,w,letter in [(ax[0],1,'a'),(ax[1],3,'b')]:
    vals=np.array([[get(rc,w,dd,c).wall_speedup for dd in ds] for c in cases]);im=a.imshow(vals,aspect='auto',vmin=0,vmax=32,cmap='Blues')
    for i in range(8):
        for j in range(3):a.text(j,i,f'{vals[i,j]:.1f}×',ha='center',va='center',fontsize=7,color='white' if vals[i,j]>20 else INK)
    a.set(xticks=range(3),xticklabels=dl,yticks=range(8),yticklabels=names);a.tick_params(length=0);title(a,letter,'Serial speedup' if w==1 else 'Three primary workers speedup')
    cb=fig.colorbar(im,ax=a,fraction=.04,pad=.03);cb.set_label('MSA wall time / method wall time',fontsize=6.2);cb.ax.tick_params(labelsize=6)
vals=np.array([[mem[(mem.dataset==dd)&(mem.case==c)].peak_rss_gb.iloc[0] for dd in ds] for c in cases]);im=ax[2].imshow(vals,aspect='auto',cmap='Blues',vmin=0,vmax=14)
for i in range(8):
    for j in range(3):ax[2].text(j,i,f'{vals[i,j]:.2f}',ha='center',va='center',fontsize=7,color='white' if vals[i,j]>8 else INK)
ax[2].set(xticks=range(3),xticklabels=dl,yticks=range(8),yticklabels=names);ax[2].tick_params(length=0);cb=fig.colorbar(im,ax=ax[2],fraction=.04,pad=.03);cb.set_label('Native graph peak RSS (GB)',fontsize=6.2);cb.ax.tick_params(labelsize=6);title(ax[2],'c','Peak graph memory with eight workers')
offsets={'F1':(6,-9),'F3':(-11,-3),'O1':(5,5),'O3':(5,5),'H1':(-4,7),'H3':(-3,-11)}
for dd,prefix,col,lab in zip(ds,['F','O','H'],[BLUE,ORANGE,GREEN],['Fountain','Organick','HT176']):
    for w,mk in [(1,'o'),(3,'s')]:
        row=get(rc,w,dd,'dbgps2');tag=prefix+str(w);ax[3].scatter(row.wall_speedup,row.cpu_ratio_vs_msa,color=col,marker=mk,s=28,label=lab if w==1 else None)
        ax[3].annotate(tag,(row.wall_speedup,row.cpu_ratio_vs_msa),xytext=offsets[tag],textcoords='offset points',fontsize=6)
ax[3].plot([0,35],[0,35],color=GRAY,ls=':',lw=.8);ax[3].set(xlim=(0,35),ylim=(0,35),xlabel='MSA wall time / DBGPS2 wall time',ylabel='MSA CPU time / DBGPS2 CPU time');ax[3].legend(loc='lower right',fontsize=6,frameon=False);title(ax[3],'d','Wall and CPU time gains differ')
fig.subplots_adjust(left=.21,right=.98,wspace=.93,hspace=.45,top=.94,bottom=.10)
save(fig,8,'speed_and_bloom')

print('SEVEN_MAIN_FIGURES_COMPLETE',flush=True)
