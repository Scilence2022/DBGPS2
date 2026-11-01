"""Six-panel speed and primary-method wall/CPU figures for new controls."""
from pathlib import Path
import sys,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R.parent/'speed_extension'));from run_bloom_matrix import CASES
DS=['fountain','nbt_file16','song_ht176'];LABELS=['DNA Fountain','Organick File 16','HT176'];COLORS=['#737373','#E69F00','#0072B2','#56B4E9','#8063B0','#B695CF','#258F6C','#72B29A','#005C51']
def main():
 s=json.loads((R/'resource_summary.json').read_text());assert len(s)==54
 get=lambda w,d,c:next(z for z in s if z['workers']==w and z['dataset']==d and z['case']==c)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axs=plt.subplots(3,2,figsize=(10.5,8.8));fig.subplots_adjust(left=.135,right=.98,bottom=.075,top=.91,hspace=.4,wspace=.6)
 for i,(d,label) in enumerate(zip(DS,LABELS)):
  for j,w in enumerate([1,3]):
   ax=axs[i,j];z=[get(w,d,c) for c in CASES];x=np.arange(9);ax.barh(x,[v['wall_mean'] for v in z],color=COLORS,xerr=[v['wall_sd'] for v in z],capsize=2,height=.66)
   ticks=['MSA','DBG rule',f'No BF, c={get(w,d,"dbgps2")["cutoff"]}',f'No BF, c={get(w,d,"dbgps2_alt_c")["cutoff"]}','B8','B16','F32','F34','B12/F34'];ax.set_yticks(x,ticks,fontsize=11.5);ax.invert_yaxis();ax.set_xscale('log');ax.set_xlim(min(v['wall_mean'] for v in z)*.65,max(v['wall_mean'] for v in z)*2.6);ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True);ax.tick_params(axis='x',labelsize=11.5);ax.set_title(f'{chr(97+i*2+j)}  {label}',loc='left',fontweight='bold',fontsize=13)
   ax.xaxis.set_major_locator(LogLocator(base=10,subs=(1,),numticks=6));ax.xaxis.set_minor_formatter(NullFormatter())
   if i==2:ax.set_xlabel('Reconstruction time (s)',fontsize=12)
   for k,v in enumerate(z):ax.annotate(f'{v["wall_speedup"]:.1f}×',(v['wall_mean']+v['wall_sd'],k),xytext=(4,0),textcoords='offset points',fontsize=11.5,ha='left',va='center')
 fig.text(.30,.982,'One serial worker',ha='center',va='top',fontweight='bold',fontsize=14);fig.text(.80,.982,'Three primary workers',ha='center',va='top',fontweight='bold',fontsize=14)
 for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig12_resource_control_speed.{ext}',dpi=300)
 plt.close(fig)
 # End-to-end wall/CPU for primary methods; all plotted quantities measured anew.
 fig,axs=plt.subplots(1,3,figsize=(10.5,4.6));fig.subplots_adjust(left=.08,right=.985,bottom=.25,top=.88,wspace=.25)
 methods=['cluster_msa_vote','dbgps_rule','dbgps2'];names=['MSA 1','DBG 1','DB2 1','MSA 3','DBG 3','DB2 3'];x=np.arange(6)
 for d,label,ax in zip(DS,LABELS,axs):
  z=[get(w,d,c) for w in [1,3] for c in methods];ax.bar(x-.2,[v['wall_mean'] for v in z],width=.38,color='#0072B2',label='Wall time',yerr=[v['wall_sd'] for v in z],capsize=2);ax.bar(x+.2,[v['cpu_mean'] for v in z],width=.38,color='#E69F00',label='User + system CPU',yerr=[v['cpu_sd'] for v in z],capsize=2);ax.set_yscale('log');ax.set_xticks(x,names,rotation=45,ha='right',fontsize=12);ax.set_title(label,fontweight='bold',fontsize=13);ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True);ax.tick_params(axis='y',labelsize=11.5)
 axs[0].set_ylabel('Time (s)');handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1.0),ncol=2,frameon=False,fontsize=12)
 for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig13_resource_control_wall_cpu.{ext}',dpi=300)
 plt.close(fig);print('RESOURCE_FIGURES_COMPLETE',flush=True)
if __name__=='__main__':main()
