"""Aggregate three repeats and generate manuscript figure/source tables."""
import csv
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
METHODS = ['cluster_msa_vote','dbgps_rule','dbgps2']
LABELS = ['Clustering + MSA + vote','DBGPS rule','DBGPS2']
COLORS = ['#9299A0','#E29545','#236F9B']


def main():
    records = json.loads((ROOT/'results.json').read_text())
    assert len(records)==18
    summary=[]
    for dataset in ['file1','file2']:
        for method in METHODS:
            rows=[r for r in records if r['dataset']==dataset and r['method']==method]
            assert len(rows)==3
            assert len({r['output_sha256'] for r in rows})==1
            times=[r['wall_seconds'] for r in rows]
            item=dict(dataset=dataset,method=method,n=3,
                      mean_seconds=statistics.mean(times),sd_seconds=statistics.stdev(times),
                      min_seconds=min(times),max_seconds=max(times),
                      correct=rows[0]['correct'],incorrect=rows[0]['incorrect'],
                      exact_recovery=rows[0]['exact_recovery'],k=rows[0]['k'],cutoff=rows[0]['cutoff'])
            for key in rows[0]['stages']:
                values=[r['stages'][key] for r in rows]
                item[key+'_mean']=statistics.mean(values)
                item[key+'_sd']=statistics.stdev(values)
            summary.append(item)
    for item in summary:
        msa=next(x for x in summary if x['dataset']==item['dataset'] and x['method']=='cluster_msa_vote')
        item['speedup_vs_msa']=msa['mean_seconds']/item['mean_seconds']
    (ROOT/'summary.json').write_text(json.dumps(summary,indent=2))
    columns=['dataset','method','n','mean_seconds','sd_seconds','min_seconds','max_seconds',
             'correct','incorrect','exact_recovery','k','cutoff','speedup_vs_msa']
    with (ROOT/'Table5_reconstruction_speed.csv').open('w') as handle:
        writer=csv.DictWriter(handle,columns,extrasaction='ignore');writer.writeheader();writer.writerows(summary)
    columns=sorted({k for x in summary for k in x})
    with (ROOT/'stage_timings.csv').open('w') as handle:
        writer=csv.DictWriter(handle,columns);writer.writeheader();writer.writerows(summary)
    flattened=[{**{k:v for k,v in r.items() if k not in ['command','stages']},**r['stages']} for r in records]
    with (ROOT/'all_repeats.csv').open('w') as handle:
        writer=csv.DictWriter(handle,sorted({k for r in flattened for k in r}));writer.writeheader();writer.writerows(flattened)
    report=['# Reconstruction speed comparison','',
            'Mean ± sample standard deviation of three serial computational repeats. '
            'Inputs are the same processed gzipped FASTA for every method. '
            'Exact recovery and wrong counts refer to unique strict strand outputs, before outer decoding.','',
            '| Input | Method | Time (s) | Exact recovery (%) | Wrong outputs | Speedup vs MSA |',
            '|---|---|---:|---:|---:|---:|']
    for row in summary:
        label=LABELS[METHODS.index(row['method'])]
        report.append(f'| {row["dataset"]} | {label} | {row["mean_seconds"]:.3f} ± {row["sd_seconds"]:.3f} | {100*row["exact_recovery"]:.4f} | {row["incorrect"]} | {row["speedup_vs_msa"]:.3f}× |')
    report.extend(['','The DBGPS comparator is its absolute-cutoff rule emulated in the same format adapter, '
                   'with a pilot-selected cutoff. DBGPS2 uses the existing c=1, r=2.55 settings. '
                   'The prefix16/MUSCLE 3.8.31 baseline is not the original full LSH pipeline. '
                   'These are implementation-specific operating-point comparisons; recovery levels differ. '
                   'The baseline invokes MUSCLE for each cluster/subset, including launch and communication costs.','',
                   'See README.md for the complete protocol, exclusions, configured worker counts and limitations. '
                   'See all_repeats.csv, stage_timings.csv, raw run logs and output SHA256 digests for the evidence.'])
    (ROOT/'BENCHMARK_REPORT.md').write_text('\n'.join(report)+'\n')
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(7.4,3.65))
    x=np.arange(2);width=.23
    for j,(method,label,color) in enumerate(zip(METHODS,LABELS,COLORS)):
        rows=[next(r for r in summary if r['dataset']==d and r['method']==method) for d in ['file1','file2']]
        loc=x+(j-1)*width
        heights=[r['mean_seconds'] for r in rows]
        axes[0].bar(loc,heights,width*.92,color=color,label=label,
                    yerr=[r['sd_seconds'] for r in rows],capsize=2,
                    error_kw={'elinewidth':.8})
        for xx,y,row in zip(loc,heights,rows):
            axes[0].text(xx,(y+row['sd_seconds'])*1.14,f'{y:.1f}',ha='center',va='bottom',fontsize=8)
        bars=axes[1].bar(loc,[100*r['exact_recovery'] for r in rows],width*.92,color=color)
        for bar,r in zip(bars,rows):
            axes[1].text(bar.get_x()+bar.get_width()/2,bar.get_height()+1.5,
                         f'{100*r["exact_recovery"]:.2f}',ha='center',fontsize=8)
    axes[0].set_yscale('log');axes[0].set_ylabel('Reconstruction time (s, log scale)')
    axes[0].set_ylim(max(.1,min(r['mean_seconds'] for r in summary)*.45),max(r['mean_seconds']+r['sd_seconds'] for r in summary)*2.5)
    axes[1].set_ylabel('Exact strand recovery (%)');axes[1].set_ylim(0,110)
    axes[1].set_yticks([0,25,50,75,100])
    for ax in axes:
        ax.set_xticks(x,['File 1\n14.30 million reads','File 2\n33.40 million reads'])
        ax.tick_params(axis='both',length=3,width=.7)
        ax.grid(axis='y',alpha=.15,linewidth=.5);ax.set_axisbelow(True)
    axes[0].text(-.18,1.04,'a',transform=axes[0].transAxes,fontsize=12,fontweight='bold')
    axes[1].text(-.18,1.04,'b',transform=axes[1].transAxes,fontsize=12,fontweight='bold')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.51,1.015),
               ncol=3,frameon=False,columnspacing=1.3,handlelength=1.2,fontsize=8.5)
    fig.subplots_adjust(left=.09,right=.985,bottom=.22,top=.84,wspace=.4)
    for extension in ['png','pdf','svg']:
        fig.savefig(ROOT/f'Fig9_reconstruction_speed.{extension}',dpi=400)
    plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
