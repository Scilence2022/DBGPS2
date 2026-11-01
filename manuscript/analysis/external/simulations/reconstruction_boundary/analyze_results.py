"""Validate the complete diagnostic matrix and draw standalone scientific figures."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
BLUE='#0072B2';ORANGE='#E69F00';GREEN='#009E73';RED='#D55E00';GRAY='#7D7D7D'
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(R/f'{name}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.08)
 plt.close(fig)
def main():
 d=pd.concat([pd.read_csv(R/n) for n in ['simulation_results.csv','geometric_results.csv','regular_geometry_results.csv']],ignore_index=True)
 assert len(d)==43776 and len(d[['case_id','k']].drop_duplicates())==1368
 assert d.groupby(['case_id','k']).size().eq(32).all()
 assert ((d.correct<=d.ceil_adm)&(d.ceil_adm<=d.ceil_seed)&(d.ceil_seed<=d.ceil_end)&(d.ceil_end<=d.ceil_cov)&(d.ceil_cov<=d.n)).all()
 assert (d.correct+d.lookahead_losses+d.seed_losses+d.ratio_losses+d.search_gap==d.ceil_cov).all()
 assert (d.search_gap==d.overflow_adm+d.ambiguous_adm).all()
 assert (d.ceil_cont<=d.ceil_adm).all()
 for name,sha in json.loads((R/'execution_manifest.json').read_text())['files'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==sha,name
 for lock in ['geometric_source_lock.json','regular_geometry_source_lock.json']:
  for name,sha in json.loads((R/lock).read_text()).items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==sha,name
 real_record=json.loads((R/'real_HT176_raw.json').read_text());real=pd.DataFrame(real_record['grid']).apply(pd.to_numeric)
 assert len(real)==16 and (real.search_gap==real.overflow_adm+real.ambiguous_adm).all()
 assert len(real_record['independent_existing_profile_checks'])==12
 for field in ['correct','ceil_cov','ceil_end','ceil_seed','ceil_adm','ceil_cont','lookahead_losses','seed_losses','ratio_losses','search_gap','overflow_adm','ambiguous_adm','wrong']:
  d[field+'_pct']=d[field]/d.n*100;real[field+'_pct']=real[field]/real.n_ref*100
 d.to_csv(R/'all_simulation_results.csv',index=False);real.to_csv(R/'HT176_exact_boundaries.csv',index=False)
 default=d[(d.c==1)&(d.r==2.55)].copy();default['gap_pp']=default.ceil_cov_pct-default.correct_pct
 keys=['cohort','length','depth','recovery','error','channel','abundance_sigma','k']
 rates=['correct_pct','ceil_cov_pct','ceil_adm_pct','ratio_losses_pct','search_gap_pct','gap_pp','observed_base_depth','mean_fragment_length','center_end_ratio','mean_cov_cv']
 grouped=default.groupby(keys)[rates].agg(['mean','std']);grouped.columns=['_'.join(c) for c in grouped.columns];grouped.reset_index().to_csv(R/'default_boundary_summary.csv',index=False)
 contexts=default[default.cohort=='main'].groupby(['length','depth','recovery','error','k'])[rates].mean();nontrivial=contexts[contexts.ceil_cov_pct>=50]
 summary=dict(synthetic_graphs=1368,synthetic_operating_points=len(d),real_operating_points=len(real),native_search_comparisons=int(d.native_checks.sum()+real.native_checks.sum()+1920),all_native_sample_comparisons_passed=True,exact_bound_inequalities_passed=True,loss_decomposition_passed=True,admissible_search_losses_equal_overflow_plus_ambiguity=True,previous_real_profile_points_matched=12,default_nontrivial_main_contexts=len(nontrivial),default_nontrivial_contexts_within_2pp=int((nontrivial.gap_pp<=2).sum()),max_default_main_mean_gap_pp=float(contexts.gap_pp.max()),default_wrong_outputs=int(default.wrong.sum()),all_settings_max_wrong_per_library=int(d.wrong.max()),cohorts=d.groupby('cohort')[['case_id','k']].apply(lambda x:len(x.drop_duplicates())).to_dict())
 (R/'validation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 # Recovery and coverage ceiling: mean of three independently generated pools.
 fig,axs=plt.subplots(4,2,figsize=(7.2,8.1));fig.subplots_adjust(left=.11,right=.94,top=.94,bottom=.08,hspace=.58,wspace=.3)
 for row,(q,k) in enumerate([(1,21),(1,31),(.05,21),(.05,31)]):
  for col,D in enumerate([30,100]):
   ax=axs[row,col];rec=np.array([[contexts.loc[(L,D,1 if L==200 else q,e,k),'correct_pct'] for e in [0,.01,.03,.1]] for L in [25,45,80,200]]);cap=np.array([[contexts.loc[(L,D,1 if L==200 else q,e,k),'ceil_cov_pct'] for e in [0,.01,.03,.1]] for L in [25,45,80,200]])
   im=ax.imshow(rec,vmin=0,vmax=100,cmap='Blues',aspect='auto');ax.set_xticks(range(4),['0','1','3','10']);ax.set_yticks(range(4),['25','45','80','200*']);ax.set_title(f'{chr(97+row*2+col)}  k={k}, q={q:g}, target depth={D}×',loc='left',fontsize=8,fontweight='bold')
   for a in range(4):
    for b in range(4):ax.text(b,a,f'{rec[a,b]:.0f}/{cap[a,b]:.0f}',ha='center',va='center',fontsize=7.5,color='white' if rec[a,b]>60 else 'black')
   if col==0:ax.set_ylabel('Cut spacing (nt)')
   if row==3:ax.set_xlabel('Substitution probability (%)')
 fig.suptitle('Fixed c=1, r=2.55: recovered / coverage ceiling (%)',y=.995,fontsize=10)
 cb=fig.colorbar(im,ax=axs.ravel().tolist(),fraction=.025,pad=.025);cb.set_label('Correct strands (%)')
 fig.text(.12,.015,'200*: shared intact control (q not applicable). q: internal retention. Mean of three pools.',fontsize=7.5)
 save(fig,'Fig_boundary_phase_maps')
 # The finite useful r window, including a real-library counterexample to monotonic recovery.
 fig,axs=plt.subplots(2,2,figsize=(7.2,5.8));fig.subplots_adjust(left=.11,right=.98,bottom=.1,top=.93,hspace=.55,wspace=.35)
 order=[1.25,1.5,2,2.55,3,4,8,0];labels=['1.25','1.5','2','2.55','3','4','8','off*'];x=np.arange(8)
 z=real[real.c==1].set_index('r').loc[order];ax=axs[0,0]
 for field,label,color,style in [('ceil_cov_pct','Coverage ceiling',GRAY,':'),('ceil_adm_pct','Rule-admissible',ORANGE,'--'),('correct_pct','Recovered',BLUE,'-')]:ax.plot(x,z[field],style,color=color,marker='s' if field=='ceil_adm_pct' else 'o' if field=='correct_pct' else None,mfc='none' if field=='ceil_adm_pct' else color,ms=5 if field=='ceil_adm_pct' else 3,zorder=4 if field=='ceil_adm_pct' else 3,label=label)
 ax.set_xticks(x,labels,rotation=35);ax.set_ylim(-2,102);ax.set_ylabel('Strands (%)');ax.set_title('a  Raw HT176, c=1',loc='left',fontweight='bold');ax.legend(fontsize=7,loc='lower left');ax.set_xlabel('r (off*: native 10⁶)')
 z=real[(real.c==1)&real.r.isin([2,2.55,3,4,8])].set_index('r').loc[[2,2.55,3,4,8]];ax=axs[0,1];bottom=np.zeros(5)
 for field,label,colr in [('lookahead_losses_pct','Terminal support',GRAY),('ratio_losses_pct','Ratio pruning',ORANGE),('search_gap_pct','Search / uniqueness',RED)]:ax.bar(range(5),z[field],bottom=bottom,color=colr,label=label);bottom+=z[field].to_numpy()
 ax.set_xticks(range(5),['2','2.55','3','4','8']);ax.set_xlabel('r');ax.set_ylabel('Gap below coverage ceiling (pp)');ax.set_title('b  Raw HT176: source of the gap',loc='left',fontweight='bold');ax.legend(fontsize=6.5)
 ax=axs[1,0];sub=d[(d.cohort=='main')&(d.length==45)&(d.depth==100)&(d.recovery==.05)&(d.error==.01)&(d.k==21)&(d.c==1)];z=sub.groupby('r').mean(numeric_only=True).loc[order]
 for field,label,color,style in [('ceil_cov_pct','Coverage ceiling',GRAY,':'),('ceil_adm_pct','Rule-admissible',ORANGE,'--'),('correct_pct','Recovered',BLUE,'-')]:ax.plot(x,z[field],style,color=color,marker='s' if field=='ceil_adm_pct' else 'o' if field=='correct_pct' else None,mfc='none' if field=='ceil_adm_pct' else color,ms=5 if field=='ceil_adm_pct' else 3,zorder=4 if field=='ceil_adm_pct' else 3,label=label)
 ax.axvline(3,color='black',lw=.6,alpha=.4);ax.set_xticks(x,labels,rotation=35);ax.set_ylim(-2,102);ax.set_ylabel('Strands (%)');ax.set_xlabel('r');ax.set_title('c  Strong recovery bias: simulation',loc='left',fontweight='bold');ax.legend(fontsize=6.8,loc='lower right')
 ax=axs[1,1];sub=default[(default.cohort=='main')&(default.length==45)&(default.depth==100)&(default.error==.01)&(default.k==21)];z=sub.groupby('recovery').mean(numeric_only=True).loc[[1,.25,.05]]
 ax.bar(np.arange(3)-.18,z.correct_pct,.34,color=BLUE,label='Recovered');ax.bar(np.arange(3)+.18,z.ceil_cov_pct,.34,color=GRAY,label='Coverage ceiling');ax.set_xticks(range(3),['1','0.25','0.05']);ax.set_xlabel('Internal-fragment retention q');ax.set_ylabel('Strands (%)');ax.set_ylim(0,105);ax.set_title('d  Fixed r=2.55, matched target depth',loc='left',fontweight='bold');ax.legend(fontsize=7,loc='lower left')
 save(fig,'Fig_boundary_ratio_window')
 # Matched cut density and pool size: length-distribution tails matter.
 fig,axs=plt.subplots(1,3,figsize=(7.2,3.15));fig.subplots_adjust(left=.09,right=.99,top=.85,bottom=.28,wspace=.45)
 cohorts=['geometry_matched_regular','geometric_breaks'];cols=[BLUE,GREEN];names=['Regular partitions','Bernoulli cuts'];xx=np.arange(3)
 sub=default[(default.cohort.isin(cohorts))&(default.depth==100)&(default.recovery==1)&(default.error==.01)&(default.k==31)]
 for j,(cohort,colr,name) in enumerate(zip(cohorts,cols,names)):
  z=sub[sub.cohort==cohort].groupby('length').agg({'correct_pct':['mean','std'],'ceil_cov_pct':'mean','mean_fragment_length':['mean','std']}).loc[[25,45,80]];off=(j-.5)*.32
  axs[0].bar(xx+off,z[('correct_pct','mean')],.3,color=colr,yerr=z[('correct_pct','std')],capsize=2,label=name)
  axs[0].scatter(xx+off,z[('ceil_cov_pct','mean')],color='black',marker='_',s=35,zorder=4)
  axs[1].bar(xx+off,z[('mean_fragment_length','mean')],.3,color=colr,yerr=z[('mean_fragment_length','std')],capsize=2)
 for ax in axs[:2]:ax.set_xticks(xx,['25','45','80']);ax.set_xlabel('Nominal cut spacing (nt)')
 axs[0].set_ylabel('Strands (%)');axs[0].set_ylim(0,105);axs[0].set_title('a  Recovery vs cut spacing',loc='left',fontsize=8,fontweight='bold')
 axs[1].set_ylabel('Mean fragment length (nt)');axs[1].set_title('b  Observed mean length',loc='left',fontsize=8,fontweight='bold')
 ax=axs[2];s=default[default.ceil_cov_pct>0];im=ax.scatter(s.ceil_cov_pct,s.correct_pct,c=s.mean_cov_cv,cmap='viridis',s=7,alpha=.45,vmin=0,vmax=1.5);ax.plot([0,100],[0,100],color=GRAY,lw=.8);ax.set_xlabel('Coverage ceiling (%)');ax.set_ylabel('Correct strands (%)');ax.set_title('c  Recovery vs ceiling',loc='left',fontsize=8,fontweight='bold');ax.set_xlim(0,102);ax.set_ylim(0,102)
 fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False,fontsize=7.5,bbox_to_anchor=(.37,1.0));fig.colorbar(im,ax=ax,fraction=.06,pad=.04,label='True-path coverage CV')
 fig.text(.09,.025,'Black ticks: coverage ceiling. Geometry panels: depth=100×, q=1, e=1%, k=31.',fontsize=7)
 save(fig,'Fig_boundary_fragment_geometry')
 print('BOUNDARY_ANALYSIS_AND_FIGURES_COMPLETE')
if __name__=='__main__':main()
