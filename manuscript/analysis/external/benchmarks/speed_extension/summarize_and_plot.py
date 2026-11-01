import csv,json,statistics
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=Path(__file__).resolve().parent
rows=json.loads((R/'results.json').read_text());assert len(rows)==27
labels={'fountain':'DNA Fountain','nbt_file16':'Organick File 16','song_ht176':'Song HT176'}
methods=['cluster_msa_vote','dbgps_rule','dbgps2'];names={'cluster_msa_vote':'Clustering + MSA + vote','dbgps_rule':'DBGPS rule','dbgps2':'DBGPS2'}
summary=[]
for d in labels:
 grouped={m:[r for r in rows if r['dataset']==d and r['method']==m] for m in methods}
 baseline=statistics.mean(r['wall_seconds'] for r in grouped['cluster_msa_vote'])
 for m in methods:
  rs=grouped[m];assert len(rs)==3 and len({x['sha256'] for x in rs})==1
  mu=statistics.mean(r['wall_seconds'] for r in rs);s=statistics.stdev(r['wall_seconds'] for r in rs)
  summary.append(dict(dataset=d,label=labels[d],method=m,method_label=names[m],time_mean=mu,time_sd=s,speedup=baseline/mu,n=3,emitted=rs[0]['emitted'],exact=rs[0]['exact'],incorrect=rs[0]['incorrect'],recovery=rs[0]['recovery'],k=rs[0]['k'],cutoff=rs[0]['cutoff'],ratio=rs[0]['ratio'],sha256=rs[0]['sha256']))
(R/'summary.json').write_text(json.dumps(summary,indent=2))
with (R/'Table6_reconstruction_speed.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
with (R/'all_repeats.csv').open('w',newline='') as f:
 fields=['dataset','method','repeat','wall_seconds','emitted','exact','incorrect','recovery','k','cutoff','ratio','sha256'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in rows)
with (R/'stage_timings.csv').open('w',newline='') as f:
 fields=['dataset','method','repeat','wall_seconds','count_seconds','seed_seconds','search_seconds','cluster_seconds','msa_vote_seconds','validation_seconds'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for row in rows:
  values={k:row[k] for k in ['dataset','method','repeat','wall_seconds']};values.update({k:row['stages'].get(k) for k in fields[4:]});w.writerow(values)
colors=['#737373','#E69F00','#0072B2'];plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(7.2,6.6));gs=fig.add_gridspec(2,2,height_ratios=[1.12,1],hspace=.65,wspace=.48);a=fig.add_subplot(gs[0,:]);b=fig.add_subplot(gs[1,0]);c=fig.add_subplot(gs[1,1])
ds=list(labels);x=np.arange(3);width=.23
for i,m in enumerate(methods):
 s=[next(z for z in summary if z['dataset']==d and z['method']==m) for d in ds];a.bar(x+(i-1)*width,[z['time_mean'] for z in s],width,yerr=[z['time_sd'] for z in s],color=colors[i],label=names[m],capsize=2)
a.set_yscale('log');a.set_ylabel('Reconstruction time (s)');a.set_xticks(x,list(labels.values()));a.set_title('a   Reconstruction after common preprocessing',loc='left',fontweight='bold');a.legend(loc='lower left',bbox_to_anchor=(0,1.05),ncol=3,frameon=False,fontsize=8);a.grid(axis='y',alpha=.2);a.set_axisbelow(True)
for i,m in enumerate(methods):
 s=[next(z for z in summary if z['dataset']==d and z['method']==m) for d in ds[1:]];b.bar(np.arange(2)+(i-1)*width,[100*z['recovery'] for z in s],width,color=colors[i]);
b.set_xticks(np.arange(2),['Organick\nFile 16','Song\nHT176']);b.set_ylim(0,105);b.set_ylabel('Exact encoded regions (%)');b.set_title('b   Exact sequence recovery',loc='left',fontweight='bold');b.grid(axis='y',alpha=.2);b.set_axisbelow(True)
f=[next(z for z in summary if z['dataset']=='fountain' and z['method']==m) for m in methods];c.bar(np.arange(3),[z['emitted'] for z in f],color=colors,width=.65)
c.set_xticks(np.arange(3),['CL–MSA\n–vote','DBGPS\nrule','DBGPS2']);c.set_ylabel('RS-valid unique seeds');c.set_ylim(0,80000);c.axhline(72000,color='black',ls='--',lw=.8);c.text(-.4,75500,'72,000 synthesized droplets',fontsize=7)
for i,z in enumerate(f):c.text(i,z['emitted']+700,f"{z['emitted']:,}",ha='center',fontsize=8)
c.set_title('c   DNA Fountain code-valid output',loc='left',fontweight='bold');c.grid(axis='y',alpha=.2);c.set_axisbelow(True)
fig.subplots_adjust(left=.105,right=.985,bottom=.08,top=.88)
for ext in ['png','pdf','svg']:fig.savefig(R/f'Fig10_extended_reconstruction_speed.{ext}',dpi=300)
print(json.dumps(summary,indent=2))
