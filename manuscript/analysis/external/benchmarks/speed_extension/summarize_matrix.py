"""Scientific figures and source tables for the complete configuration matrix."""
import csv,json,statistics
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from run_bloom_matrix import CASES
R=Path(__file__).resolve().parent
DATASETS=['fountain','nbt_file16','song_ht176'];LABELS=['DNA Fountain','Organick File 16','HT176']
COLORS=['#737373','#E69F00','#0072B2','#56B4E9','#8063B0','#B695CF','#258F6C','#72B29A','#005C51']
NAMES=['Clustering + MSA + vote','DBGPS rule','DBGPS2 no Bloom','DBGPS2 no Bloom alternate c','DBGPS2 B8','DBGPS2 B16','DBGPS2 F32','DBGPS2 F34','DBGPS2 B12 F34']
def main():
 rows=json.loads((R/'matrix_results.json').read_text());assert len(rows)==81
 fountain_checks={x['case']:x for x in json.loads((R/'matrix_validation.json').read_text())['fountain']}
 summary=[]
 for dataset,label in zip(DATASETS,LABELS):
  baseline=statistics.mean(r['wall_seconds'] for r in rows if r['dataset']==dataset and r['case']=='cluster_msa_vote')
  for case,name in zip(CASES,NAMES):
   rs=[r for r in rows if r['dataset']==dataset and r['case']==case];assert len(rs)==3 and len({r['sha256'] for r in rs})==1
   z=rs[0];mu=statistics.mean(r['wall_seconds'] for r in rs);rss=[r['reported_graph_rss_gb'] for r in rs if r['reported_graph_rss_gb'] is not None]
   summary.append(dict(dataset=dataset,label=label,case=case,method_label=name,time_mean=mu,time_sd=statistics.stdev(r['wall_seconds'] for r in rs),speedup=baseline/mu,n=3,emitted=z['emitted'],exact=z['exact'],incorrect=z['incorrect'],recovery=z['recovery'],k=z['k'],cutoff=z['cutoff'],ratio=z['ratio'],query_bits=z['query_bits'],filter_shift=z['filter_shift'],query_bytes=z.get('query_bytes',0),peak_rss_gb=max(rss) if rss else None,sha256=z['sha256'],payload_inconsistent=fountain_checks[case].get('payload_mismatches_against_reencoded_archive') if dataset=='fountain' else None))
 (R/'matrix_summary.json').write_text(json.dumps(summary,indent=2))
 with (R/'Table6_configuration_comparison.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
 fields=['dataset','case','repeat','wall_seconds','emitted','exact','incorrect','recovery','k','cutoff','ratio','query_bits','filter_shift','query_bytes','maximum_rss_bytes','reported_graph_rss_gb','sha256']
 with (R/'matrix_all_repeats.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r.get(k) for k in fields} for r in rows)
 fields=['dataset','case','repeat','count_seconds','bloom_seconds','seed_seconds','search_seconds','cluster_seconds','msa_vote_seconds','validation_seconds']
 with (R/'matrix_stage_timings.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for r in rows:w.writerow(dict(dataset=r['dataset'],case=r['case'],repeat=r['repeat'],**{k:r['stages'].get(k,0) for k in fields[3:]}))
 get=lambda d,c:next(z for z in summary if z['dataset']==d and z['case']==c)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axes=plt.subplots(3,1,figsize=(7.2,6.6));fig.subplots_adjust(left=.12,right=.985,bottom=.08,top=.95,hspace=.72)
 for i,(d,label,ax) in enumerate(zip(DATASETS,LABELS,axes)):
  zs=[get(d,c) for c in CASES];x=np.arange(len(zs));ax.bar(x,[z['time_mean'] for z in zs],color=COLORS,width=.72,yerr=[z['time_sd'] for z in zs],capsize=2)
  ticks=['CL–MSA\n–vote','DBGPS\nrule',f'No BF\nc={get(d,"dbgps2")["cutoff"]}',f'No BF\nc={get(d,"dbgps2_alt_c")["cutoff"]}','B8','B16','F32','F34','B12 +\nF34'];ax.set_xticks(x,ticks,fontsize=8);ax.set_yscale('log');ax.set_ylabel('Reconstruction time (s)');ax.set_title(f'{chr(97+i)}   {label}',loc='left',fontweight='bold');ax.set_ylim(top=max(z['time_mean'] for z in zs)*2.5);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
  for j,z in enumerate(zs):ax.text(j,z['time_mean']*1.13,f'{z["speedup"]:.1f}×',ha='center',va='bottom',fontsize=7)
 for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig10_reconstruction_speed_Bloom_matrix.{ext}',dpi=300)
 plt.close(fig)
 # Graph-only RSS and output-yield tradeoff, relative to primary no-Bloom point.
 graphcases=CASES[1:];mem=np.array([[get(d,c)['peak_rss_gb'] for d in DATASETS] for c in graphcases]);delta=np.array([[(get(d,c)['emitted'] if d=='fountain' else get(d,c)['exact'])-(get(d,'dbgps2')['emitted'] if d=='fountain' else get(d,'dbgps2')['exact']) for d in DATASETS] for c in graphcases])
 fig,(a,b)=plt.subplots(1,2,figsize=(7.2,5.1),sharey=True);fig.subplots_adjust(left=.22,right=.97,bottom=.23,top=.91,wspace=.35)
 im=a.imshow(mem,aspect='auto',cmap='Blues');limit=max(1,float(np.max(np.abs(delta))));im2=b.imshow(delta,aspect='auto',cmap='RdBu',vmin=-limit,vmax=limit)
 yl=['DBGPS rule','No Bloom, primary c','No Bloom, alternate c','Query B8','Query B16','Counting F32','Counting F34','B12 + F34'];a.set_yticks(np.arange(8),yl,fontsize=8)
 for ax in [a,b]:ax.set_xticks(range(3),['Fountain','Organick\nFile 16','HT176'],fontsize=8);ax.tick_params(axis='both',length=0)
 a.set_title('a   Graph peak RSS (GB)',loc='left',fontweight='bold');b.set_title('b   Change in output yield',loc='left',fontweight='bold')
 for i in range(8):
  for j in range(3):
   a.text(j,i,f'{mem[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if mem[i,j]>.65*mem.max() else 'black');b.text(j,i,f'{int(delta[i,j]):+,}',ha='center',va='center',fontsize=8,color='white' if abs(delta[i,j])>.65*limit else 'black')
 cb=fig.colorbar(im,ax=a,location='bottom',fraction=.06,pad=.2);cb.ax.tick_params(labelsize=7);cb=fig.colorbar(im2,ax=b,location='bottom',fraction=.06,pad=.2);cb.ax.tick_params(labelsize=7)
 fig.text(.22,.015,'Yield change versus the primary no-Bloom point: RS-valid seeds for Fountain;\nexact encoded regions for Organick and HT176. RSS is graph-process memory.',fontsize=8)
 for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig11_Bloom_memory_and_recovery.{ext}',dpi=300)
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
