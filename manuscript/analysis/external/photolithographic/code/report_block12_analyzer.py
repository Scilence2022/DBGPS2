#!/usr/bin/env python3
"""Interpret actual DBGPS-analyzer profiles against File 3 recovery evidence."""
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/analyzer_block12'
INNER=DATA/'results/inner_rs'
UNION=INNER/'parameter_union'
N=196596
B=16383
KS=[13,15,17,18,19,20,21,25,31]
COMP=str.maketrans('ACGT','TGCA')

def load_counts(name,k,seqs):
    path=OUT/(name+'.cov_details');cache=OUT/(name+'.ordered_counts.npz')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    if cache.exists():
        with np.load(cache) as data:
            if str(data['source_sha256'])==digest:return data['coverage'].copy()
    coverage=np.empty((N,61-k),dtype=np.uint16);repeated=0
    with path.open() as handle:
        for i,line in enumerate(handle):
            field,text=line.split('\t');assert int(field)==i
            values=np.fromstring(text,sep=' ',dtype=np.uint16)
            if len(values)==61-k:coverage[i]=values
            else:
                # Standalone profiles deduplicate canonical keys in first-occurrence
                # order. Expand repeated keys back to their physical positions.
                lookup={};cursor=0
                for pos in range(61-k):
                    key=seqs[i][pos:pos+k];key=min(key,key.translate(COMP)[::-1])
                    if key not in lookup:lookup[key]=values[cursor];cursor+=1
                    coverage[i,pos]=lookup[key]
                assert cursor==len(values);repeated+=1
        assert i+1==N
    np.savez_compressed(cache,coverage=coverage,source_sha256=digest,repeated_reference_rows=repeated)
    return coverage

def canonical_keys(seqs,k):
    dna=np.frombuffer(''.join(seqs).encode(),dtype=np.uint8).reshape(N,60)
    table=np.zeros(256,dtype=np.uint64)
    for j,c in enumerate(b'ACGT'):table[c]=j
    dna=table[dna];forward=np.zeros(N,dtype=np.uint64);reverse=forward.copy()
    mask=np.uint64((1<<(2*k))-1);keys=np.empty((N,61-k),dtype=np.uint64)
    for pos in range(60):
        forward=((forward<<np.uint64(2))|dna[:,pos])&mask
        reverse=(reverse>>np.uint64(2))|((np.uint64(3)-dna[:,pos])<<np.uint64(2*(k-1)))
        if pos>=k-1:keys[:,pos-k+1]=np.minimum(forward,reverse)
    alphabet=str.maketrans('ACGT','0123')
    for i in [0,1,16382,16383,180213,196595]:
        for pos in range(61-k):
            text=seqs[i][pos:pos+k];canonical=min(text,text.translate(COMP)[::-1])
            assert int(keys[i,pos])==int(canonical.translate(alphabet),4)
    return keys

def read_tool_metrics(name):
    with (OUT/(name+'.SmKdKn')).open() as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    return [{k:None if v=='nan' else float(v) if k in ['Ratio','Sm','Kd','Kn'] else int(v) for k,v in r.items()} for r in rows]

def profile_metrics(c,unique_indices,k):
    rows=[]
    for b in range(12):
        x=c[b*B:(b+1)*B];minimum=x.min(axis=1);means=x.mean(axis=1)
        distinct=x.ravel()[unique_indices[b]]
        rows.append(dict(block=b+1,k=k,strands=B,complete_gt0=int((minimum>0).sum()),
            complete_gt1=int((minimum>1).sum()),complete_gt3=int((minimum>3).sum()),
            all_kmers_zero=int((x.max(axis=1)==0).sum()),
            median_strand_mean_coverage=float(np.median(means)),mean_window_coverage=float(x.mean()),
            strand_mean_coverage_quantiles=[float(z) for z in np.quantile(means,[0,.1,.25,.5,.75,.9,1])],
            missing_window_fraction=float((x==0).mean()),weak_window_fraction=float((x<=1).mean()),
            distinct_reference_kmers=len(distinct),kd_gt0=float((distinct==0).mean()),kd_gt1=float((distinct<=1).mean()),
            position_zero_fraction=(x==0).mean(axis=0).tolist(),
            position_weak_fraction=(x<=1).mean(axis=0).tolist(),
            position_mean_coverage=x.mean(axis=0).tolist(),position_median_coverage=np.median(x,axis=0).tolist()))
    return rows

def ratio_allowed(c,k,r):
    anchor=54-k
    prev=np.concatenate([c[:,1:anchor+1],c[:,anchor:-1]],axis=1)
    nxt=np.concatenate([c[:,:anchor],c[:,anchor+1:]],axis=1)
    return np.all((nxt>=np.floor(prev/r))&(nxt<=np.minimum(16383,np.floor(prev*r))),axis=1)

def read_overflows(k,c,r):
    path=DATA/f'results/file3.k{k}.c{c}.r{r:.2f}.status.tsv'
    with path.open() as handle:
        assert next(handle).strip()=='id\tseeds\tcomplete_paths\tvalid_paths\tstatus'
        rows=np.fromiter((line.rstrip().endswith('\toverflow') for line in handle),dtype=bool,count=N)
        assert len(rows)==N and not handle.read()
    return rows

def main():
    seqs=[s for s in (DATA/'sources/File3_synthesized.fasta').read_text().splitlines() if s and not s.startswith('>')]
    assert len(seqs)==N and all(len(s)==60 for s in seqs)
    with (UNION/'per_index.tsv').open() as handle:index_rows=list(csv.DictReader(handle,delimiter='\t'))
    recovered=np.array([bool(int(r['candidate_success'])) for r in index_rows])
    success_matrix=np.loadtxt(UNION/'candidate_success_matrix.tsv.gz',skiprows=1,usecols=range(1,180),dtype=np.uint8)
    assert success_matrix.shape==(N,179) and np.array_equal(recovered,np.any(success_matrix>0,axis=1))
    with (UNION/'parameters.tsv').open() as handle:
        grid=[dict(parameter_id=int(r['parameter_id']),k=int(r['k']),cutoff=int(r['cutoff']),ratio=float(r['ratio'])) for r in csv.DictReader(handle,delimiter='\t')]
    plan=json.loads((OUT/'plan.json').read_text());profiles={};counts={};unique={};global_unique={};design=[];validations=[]
    for run in plan:
        name=f"{run['dataset']}.k{run['k']}.L{run['limit']}";k=run['k']
        metadata=json.loads((OUT/(name+'.run.json')).read_text());assert metadata['exit_code']==0
        c=load_counts(name,k,seqs);counts[name]=c
        if k not in unique:
            keys=canonical_keys(seqs,k)
            unique[k]=[np.unique(keys[b*B:(b+1)*B].ravel(),return_index=True)[1] for b in range(12)]
            global_unique[k]=np.unique(keys.ravel(),return_index=True)[1]
        profiles[name]=dict(**run,per_block=profile_metrics(c,unique[k],k),global_tool_metrics=read_tool_metrics(name))
        distinct_global=c.ravel()[global_unique[k]]
        for tool in profiles[name]['global_tool_metrics']:
            assert tool['Paths']==int((c.min(axis=1)>tool['Coverage']).sum())
            assert tool['Exist']==int((distinct_global>tool['Coverage']).sum())
            assert tool['Lost']==int((distinct_global<=tool['Coverage']).sum())
        if run['dataset']=='processed':
            old=np.loadtxt(DATA/f'results/file3.k{k}.truth_support.tsv',skiprows=1,usecols=1,dtype=np.uint16)
            minimum=c.min(axis=1);expected=np.where(minimum>1,minimum,0)
            assert np.array_equal(old,expected),name
            validations.append(dict(name=name,all_reference_min_counts_match_original_filtered_counter=True))
        print('PROFILE',name,'block12 complete>1',profiles[name]['per_block'][11]['complete_gt1'],flush=True)
    for k in [17,19,20]:
        processed=counts[f'processed.k{k}.L60'];raw=counts[f'raw.k{k}.L200']
        assert np.all(raw>=processed)
    assert np.all(counts['raw.k19.L200']>=counts['raw.k19.L60'])
    assert np.all(counts['raw.k19.L60']>=counts['processed.k19.L60'])
    # Exact reference-path feasibility under the original directed, rounded rule.
    feasible_any=np.zeros(N,dtype=bool);survived_any=feasible_any.copy();coverage_any=feasible_any.copy()
    graph_rows=[];graph_by_k=[];best_masks=None
    for k in KS:
        c=counts[f'processed.k{k}.L60'];minimum=c.min(axis=1)
        ratios={r['ratio']:ratio_allowed(c,k,r['ratio']) for r in grid if r['k']==k}
        k_feasible=np.zeros(N,dtype=bool);k_survived=k_feasible.copy();k_coverage=k_feasible.copy()
        k_columns=[]
        for r in [r for r in grid if r['k']==k]:
            p=r['parameter_id'];eligible=minimum>r['cutoff'];feasible=eligible&ratios[r['ratio']]
            overflow=read_overflows(k,r['cutoff'],r['ratio']);survived=feasible&~overflow
            zero=success_matrix[:,p]==1
            assert np.array_equal(survived,zero),f'Exact path enumeration disagreement: {r}'
            k_coverage|=eligible;k_feasible|=feasible;k_survived|=survived;k_columns.append(p)
            graph_rows.append(dict(**r,per_block=[dict(block=b+1,coverage_eligible=int(eligible[b*B:(b+1)*B].sum()),
                ratio_eligible=int(feasible[b*B:(b+1)*B].sum()),
                feasible_but_overflow=int((feasible&overflow)[b*B:(b+1)*B].sum()),
                exact_unmodified_candidates=int(survived[b*B:(b+1)*B].sum()),
                exact_after_rs=int((success_matrix[b*B:(b+1)*B,p]>0).sum())) for b in range(12)]))
            if k==19 and r['cutoff']==1 and r['ratio']==3.5:
                best_masks=dict(minimum=minimum,eligible=eligible,feasible=feasible,overflow=overflow,survived=survived,
                                recovered=success_matrix[:,p]>0)
        feasible_any|=k_feasible;survived_any|=k_survived;coverage_any|=k_coverage
        same_k=np.any(success_matrix[:,k_columns]>0,axis=1)
        graph_by_k.append(dict(k=k,per_block=[dict(block=b+1,coverage_eligible_any=int(k_coverage[b*B:(b+1)*B].sum()),
            ratio_eligible_any=int(k_feasible[b*B:(b+1)*B].sum()),
            unmodified_candidate_any=int(k_survived[b*B:(b+1)*B].sum()),
            correct_after_rs_any=int(same_k[b*B:(b+1)*B].sum())) for b in range(12)]))
        print('GRAPH VALIDATED',k,'all parameter exact-path predictions match observed zero-error candidates',flush=True)
    assert np.array_equal(survived_any,np.any(success_matrix==1,axis=1))
    assert best_masks is not None
    best=best_masks;losses=[]
    for b in range(12):
        sl=slice(b*B,(b+1)*B);failed=~best['recovered'][sl];minimum=best['minimum'][sl]
        categories=dict(missing_kmer=int((failed&(minimum==0)).sum()),
            singleton_bottleneck=int((failed&(minimum==1)).sum()),
            rounded_ratio_rejected=int((failed&best['eligible'][sl]&~best['feasible'][sl]).sum()),
            feasible_but_overflow=int((failed&best['feasible'][sl]&best['overflow'][sl]).sum()))
        assert sum(categories.values())==int(failed.sum())
        losses.append(dict(block=b+1,failed_after_rs=int(failed.sum()),**categories))
    union_losses=[]
    for b in range(12):
        sl=slice(b*B,(b+1)*B);failed=~recovered[sl]
        categories=dict(no_coverage_eligible_setting=int((failed&~coverage_any[sl]).sum()),
            only_ratio_rejected=int((failed&coverage_any[sl]&~feasible_any[sl]).sum()),
            feasible_only_in_overflow_searches=int((failed&feasible_any[sl]&~survived_any[sl]).sum()))
        assert sum(categories.values())==int(failed.sum())
        union_losses.append(dict(block=b+1,failed_after_rs=int(failed.sum()),**categories))
    # Locate weak support by the physical regions touched by each k-mer window.
    c19=counts['processed.k19.L60'];tail=counts['raw.k19.L200'];rawcap=counts['raw.k19.L60']
    position_rows=[]
    for b in range(12):
        for pos in range(42):
            v=c19[b*B:(b+1)*B,pos];rawv=tail[b*B:(b+1)*B,pos]
            position_rows.append(dict(block=b+1,k=19,position_0based=pos,window_end_inclusive=pos+18,
                touches_payload=pos<42,touches_index=pos<=53 and pos+18>=42,touches_parity=pos+18>=54,
                processed_missing=int((v==0).sum()),processed_weak=int((v<=1).sum()),
                processed_median=float(np.median(v)),raw_missing=int((rawv==0).sum()),raw_weak=int((rawv<=1).sum()),raw_median=float(np.median(rawv))))
    for b in range(12):
        seq=seqs[b*B:(b+1)*B];gc=np.array([sum(x in 'GC' for x in s)/60 for s in seq])
        hp=np.array([max(len(list(g)) for _,g in itertools.groupby(s)) for s in seq])
        design.append(dict(block=b+1,mean_gc=float(gc.mean()),median_gc=float(np.median(gc)),
            homopolymer_ge6=int((hp>=6).sum()),homopolymer_ge8=int((hp>=8).sum()),max_homopolymer=int(hp.max())))
    links={str(k):json.loads((OUT/f'links.k{k}.json').read_text()) for k in KS}
    summary=dict(tool='DBGPS-analyzer and DBGPS-links',provenance=json.loads((OUT/'provenance.json').read_text()),
        target_block=12,index_range_0based=[180213,196595],reference_strands=N,block_size=B,
        read_sets=dict(processed=92337721,raw=195619515),profiles=profiles,
        best_single_configuration_loss_breakdown=losses,union_loss_breakdown=union_losses,
        graph_feasibility_by_k=graph_by_k,design=design,links=links,
        validation=dict(processed_counter_checks=validations,exact_reference_path_configurations=len(graph_rows),
            exact_reference_path_predictions_equal_observed_zero_error_candidates=True,
            all_13_cli_sm_and_distinct_kd_counts_match_expanded_profiles=True,
            raw_counts_dominate_processed=True,raw_full_counts_dominate_raw_capped=True),
        definitions=dict(complete_gt1='All physical reference k-mer windows have saturated pooled count >1; no ratio or ambiguity constraint',
            kd='Distinct canonical k-mers within the named block, missing from the pooled reads at the named count cutoff',
            coverage='Pooled canonical k-mer occurrence count, not assigned per-strand read depth',
            feasible='Exact unmodified reference path satisfies the adapter directed integer-rounded coverage rule',
            overflow='Original 1050-live-path search limit; the index search discards its full candidate pool on overflow',
            union_loss='Classifies absence of an unmodified exact path among failed RS-recovery indices; it does not enumerate every correctable one-symbol neighbour',
            links='Counted by DBGPS-links, sum(count-1) for keys shared by distinct reference strands; count cap1023',
            noise='Global full-read-index noise only for processed k19; Bloom runs do not estimate off-reference noise'))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'graph_parameter_diagnostics.json').write_text(json.dumps(graph_rows,indent=2))
    write_table(OUT/'all_blocks.tsv',[dict(dataset=name,**{k:v for k,v in row.items() if not isinstance(v,list)}) for name,pr in profiles.items() for row in pr['per_block']])
    write_table(OUT/'position_profiles.k19.tsv',position_rows)
    write_table(OUT/'best_parameter_loss_breakdown.tsv',losses)
    write_table(OUT/'union_loss_breakdown.tsv',union_losses)
    with (OUT/'block12_per_index.tsv').open('w') as handle:
        names=['index','recovered_in_179_parameter_union','k19_missing_windows','k19_singleton_windows','k19_mean_coverage',
               'k19_best_reference_path_eligible','k19_best_search_overflow','k19_best_correct_after_rs',
               'raw_k19_L60_missing_windows','raw_k19_L200_missing_windows']
        for k in KS:names += [f'k{k}_missing_windows',f'k{k}_min_coverage']
        writer=csv.writer(handle,delimiter='\t');writer.writerow(names)
        for i in range(11*B,N):
            row=[i,int(recovered[i]),int((c19[i]==0).sum()),int((c19[i]==1).sum()),round(float(c19[i].mean()),6),
                int(best['feasible'][i]),int(best['overflow'][i]),int(best['recovered'][i]),
                int((rawcap[i]==0).sum()),int((tail[i]==0).sum())]
            for k in KS:
                c=counts[f'processed.k{k}.L60'][i];row += [int((c==0).sum()),int(c.min())]
            writer.writerow(row)
    plot(summary)
    write_report(summary)
    files=['summary.json','REPORT.zh-CN.md','all_blocks.tsv','position_profiles.k19.tsv','block12_per_index.tsv',
           'best_parameter_loss_breakdown.tsv','union_loss_breakdown.tsv','graph_parameter_diagnostics.json','diagnostics.png','diagnostics.svg']
    (OUT/'validation.json').write_text(json.dumps(dict(**summary['validation'],sha256={f:hashlib.sha256((OUT/f).read_bytes()).hexdigest() for f in files}),indent=2))
    print('BLOCK 12 ANALYSIS COMPLETE',flush=True)

def write_table(path,rows):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)

def plot(s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    x=np.arange(1,13);p=s['profiles']['processed.k19.L60']['per_block'];raw=s['profiles']['raw.k19.L200']['per_block']
    axes[0,0].bar(x-.18,[r['complete_gt1'] for r in p],.36,label='Processed reads, k=19')
    axes[0,0].bar(x+.18,[r['complete_gt1'] for r in raw],.36,label='All raw reads, k=19')
    axes[0,0].axhline(10977,color='black',ls='--',label='Outer RS erasure threshold')
    axes[0,0].set(xlabel='Block',ylabel='Strands with all k-mers counted >=2',xticks=x,title='A. Exact-reference graph coverage');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(KS,[s['profiles'][f'processed.k{k}.L60']['per_block'][11]['complete_gt1'] for k in KS],'-o',label='Block 12: coverage only')
    axes[0,1].plot(KS,[r['per_block'][11]['unmodified_candidate_any'] for r in s['graph_feasibility_by_k']],'-o',label='Unmodified candidates: all tested c/r')
    axes[0,1].plot(KS,[r['per_block'][11]['correct_after_rs_any'] for r in s['graph_feasibility_by_k']],'-o',label='Correct after inner RS: all tested c/r')
    axes[0,1].axhline(10977,color='black',ls='--');axes[0,1].set(xlabel='k',ylabel='Block 12 strands',title='B. Coverage versus actual search');axes[0,1].legend(fontsize=8)
    a=np.array([r['position_weak_fraction'] for r in p])
    im=axes[1,0].imshow(a,aspect='auto',origin='lower',extent=[-.5,41.5,.5,12.5],vmin=0,vmax=max(.5,a.max()),cmap='magma')
    axes[1,0].set(xlabel='k-mer start position (0-based), k=19',ylabel='Block',yticks=x,title='C. Fraction of windows with count <=1');fig.colorbar(im,ax=axes[1,0],label='Fraction')
    losses=s['best_single_configuration_loss_breakdown'];bottom=np.zeros(12)
    for key,label in [('missing_kmer','Missing k-mer'),('singleton_bottleneck','Singleton bottleneck'),('rounded_ratio_rejected','Coverage-ratio rejection'),('feasible_but_overflow','Path-limit overflow')]:
        vals=np.array([r[key] for r in losses]);axes[1,1].bar(x,vals,bottom=bottom,label=label);bottom+=vals
    axes[1,1].set(xlabel='Block',ylabel='Unrecovered reference strands',xticks=x,title='D. Losses at k19 / c1 / r3.5');axes[1,1].legend(fontsize=8)
    fig.savefig(OUT/'diagnostics.png',dpi=180);fig.savefig(OUT/'diagnostics.svg');plt.close(fig)

def write_report(s):
    profiles=s['profiles'];p=profiles['processed.k19.L60']['per_block'];raw=profiles['raw.k19.L200']['per_block'];cap=profiles['raw.k19.L60']['per_block']
    b=p[11];r=raw[11];c=cap[11];loss=s['best_single_configuration_loss_breakdown'][11]
    lines=['**File 3 第 12 块：DBGPS-analyzer 诊断与其余 11 块对照**','',
        f"第 12 块为零基索引 **180213–196595**，共 16383 条。使用用户指定仓库的实际 `DBGPS-analyzer` 和 `DBGPS-links`，从当前源码独立编译，未修改工具源码。所有 12 块共同计数并分别汇总。目标保留完整 60 nt，不裁剪引物。",'',
        f"核心发现：处理后 reads、k=19 时，第 12 块所有 k-mer 都至少出现两次的目标链只有 **{b['complete_gt1']:,}** 条；逐链平均 k-mer 覆盖的中位数为 **{b['median_strand_mean_coverage']:.2f}**。完整原始 reads 提升至 **{r['complete_gt1']:,}** 条和 **{r['median_strand_mean_coverage']:.2f}**，但该 k 下仍低于外码只含擦除时所需的 10977 条。",'',
        '覆盖指全部混合 reads 对目标 canonical k-mer 的计数，包含跨链共享和重复出现贡献，不等于已经归属到每条目标链的 read 深度。有效正确 k-mer 支持低这一结论有直接证据；单凭这些统计，尚不能确定原始分子丰度低与测序/合成错误较多各占多少。','',
        '**所有 12 块的 k=19 覆盖对照**','',
        '| 块 | 处理后全部 k-mer >0 | 处理后全部 k-mer >1 | 原始全部 k-mer >1 | 处理后覆盖中位数 | 原始覆盖中位数 |',
        '|---:|---:|---:|---:|---:|---:|']
    for a,z in zip(p,raw):lines.append(f"| {a['block']} | {a['complete_gt0']:,} | {a['complete_gt1']:,} | {z['complete_gt1']:,} | {a['median_strand_mean_coverage']:.2f} | {z['median_strand_mean_coverage']:.2f} |")
    lines+=['','表中完整覆盖是参考已知条件下的图支持诊断，尚未施加覆盖比、搜索路径上限或候选唯一性要求，不是已实现的盲恢复数。','',
        '**第 12 块：长度筛选与 60 nt 截短的影响**','',
        '| 输入 | reads 数 | 所有 k-mer >0 | 所有 k-mer >1 | 覆盖中位数 |', '|---|---:|---:|---:|---:|']
    for label,reads,row in [('长度55–70，再截至60nt',92337721,b),('全部原始reads，最多取60nt',195619515,c),('全部原始reads，保留实际长度',195619515,r)]:
        lines.append(f"| {label} | {reads:,} | {row['complete_gt0']:,} | {row['complete_gt1']:,} | {row['median_strand_mean_coverage']:.2f} |")
    lines += ['',f"在 k=19、计数>1、暂不限制覆盖比时：取消长度筛选新增完整支持 **{c['complete_gt1']-b['complete_gt1']:,}** 条；进一步保留第 60 nt 之后的碱基新增 **{r['complete_gt1']-c['complete_gt1']:,}** 条。后者也可包含从较后位置开始的目标 k-mer 支持，不能仅解释为末端校验区被截掉。",'',
        f"第12块逐位置零覆盖窗口占比为 **{100*min(b['position_zero_fraction']):.2f}%–{100*max(b['position_zero_fraction']):.2f}%**，起始窗口为 {100*b['position_zero_fraction'][0]:.2f}%，末端窗口为 {100*b['position_zero_fraction'][-1]:.2f}%。低支持跨越 payload、index、parity 相关窗口；一个缺失窗口不能直接定位到某个错误碱基。",'',
        '**第 12 块：k 的变化与真实搜索之间的差距**','',
        '| k | 完整覆盖：计数>0 | 完整覆盖：计数>1 | 各 c/r 中至少有一组允许真值路径 | 真值路径且未溢出的并集 | RS 后正确候选并集 |',
        '|---:|---:|---:|---:|---:|---:|']
    for row in s['graph_feasibility_by_k']:
        k=row['k'];a=profiles[f'processed.k{k}.L60']['per_block'][11];g=row['per_block'][11]
        lines.append(f"| {k} | {a['complete_gt0']:,} | {a['complete_gt1']:,} | {g['ratio_eligible_any']:,} | {g['unmodified_candidate_any']:,} | {g['correct_after_rs_any']:,} |")
    lines+=['','各行使用此前实际测试的该 k 对应 c/r 网格，不是相同的完整笛卡尔网格。小 k 的高覆盖可能来自跨链共享及噪声碰撞；它不能保证可唯一组装。真值路径允许但未输出的情况，结合原始搜索的 overflow 标记检查。RS 可从有一个符号错误的另一条路径修回真值，因此纠错后正确候选数可以高于“未经修正真值路径”的数量。','',
        '**最佳单参数 k=19、c=1、r=3.5 下的未恢复原因**','',
        f"该配置第 12 块 RS 后未恢复 {loss['failed_after_rs']:,} 条，按互斥条件分为：",'',
        f"- 至少一个真值 k-mer 计数为 0：**{loss['missing_kmer']:,}** 条。",
        f"- 没有零计数，但最弱 k-mer 只有 1 次，被 c=1 排除：**{loss['singleton_bottleneck']:,}** 条。",
        f"- 所有 k-mer 均超过 c，但真值路径被整数取整后的覆盖比规则拒绝：**{loss['rounded_ratio_rejected']:,}** 条。",
        f"- 真值完整路径满足计数与覆盖比规则，但该索引搜索超过 1050 条活跃路径、整组被丢弃：**{loss['feasible_but_overflow']:,}** 条。",'',
        f"两类计数不足共 {loss['missing_kmer']+loss['singleton_bottleneck']:,} 条，占该配置失败数的 {100*(loss['missing_kmer']+loss['singleton_bottleneck'])/loss['failed_after_rs']:.2f}%。第6块则有 {s['best_single_configuration_loss_breakdown'][5]['feasible_but_overflow']:,} 条合格真值路径因溢出丢失，说明不同块的主要限制并不相同。",'',
        '以上从该参数 RS 后仍失败的索引中计数，解释未经修正真值路径在哪一步受阻；不是把 RS 可纠正邻居的所有可能路径逐一穷举。四项之和严格等于该配置的未恢复数。所有 12 块的同类分解另存为表。','',
        '**跨全部 179 组参数仍失败的第 12 块索引**','']
    u=s['union_loss_breakdown'][11]
    lines += [f"此前正确候选并集为 6927 条，仍失败 {u['failed_after_rs']:,} 条。按未经修正真值路径在全部网格中的可达性分类：没有任何计数合格参数 {u['no_coverage_eligible_setting']:,} 条；仅能满足计数、始终不满足覆盖比 {u['only_ratio_rejected']:,} 条；存在完整合格路径但所有相关搜索均溢出 {u['feasible_only_in_overflow_searches']:,} 条。",'',
        '该跨网格分类包含很小的 k；当短 k-mer 被噪声或其他链充分覆盖时，即使真实分子证据不足，参考序列也可能在图上形成路径。因此“因溢出而丢弃合格真值路径”不证明放大路径上限就能可靠恢复这些链。还需评价候选数量、RS 误纠正和唯一选择。','',
        '**原始 reads 的低 k 对照**','', '| k | 处理后完整覆盖 >1 | 原始完整覆盖 >1 | 原始完整覆盖 >0 |','|---:|---:|---:|---:|']
    for k in [17,19,20]:
        a=profiles[f'processed.k{k}.L60']['per_block'][11];z=profiles[f'raw.k{k}.L200']['per_block'][11]
        lines.append(f"| {k} | {a['complete_gt1']:,} | {z['complete_gt1']:,} | {z['complete_gt0']:,} |")
    noise=profiles['processed.k19.L60']['global_tool_metrics'][1]
    lines += ['',f"完整索引的全库 k=19、c=1 统计：目标存在 k-mer **{noise['Exist']:,}** 个，非目标噪声 k-mer **{noise['Noise']:,}** 个，`Kn={noise['Kn']:.6f}`。这是不同 k-mer 的数量比，不是碱基错误率；不能把其他块的正常链相对某一单块称为噪声，因此此处只报告全库 Kn。其余计数使用工具的 exact-reference Bloom 后端，目标覆盖保持精确，未估计全局 Kn。",'',
        '**参考设计的链间共享**','', '| k | 全库交叉链接 | 第12块内部交叉链接 | 其他单块内部范围 |','|---:|---:|---:|---|']
    for k in KS:
        a=s['links'][str(k)];v=[x['cross_links'] for x in a[1:12]]
        lines.append(f"| {k} | {a[0]['cross_links']:,} | {a[12]['cross_links']:,} | {min(v):,}–{max(v):,} |")
    d=s['design'][11];others=s['design'][:11]
    lines+=['','链接由 `DBGPS-links -p 0` 计算，为被多个目标链共享的 canonical k-mer 的额外链成员数之和。全库与单块大小不同，不能直接按总数比较；单块之间包含相同的 16383 条链。工具单键计数上限为 1023。此项只测设计共享，不测含错 reads 形成的额外图分支。','',
        f"参考设计方面，第12块平均 GC 为 {100*d['mean_gc']:.2f}%，含至少8 nt同聚物的链有 {d['homopolymer_ge8']} 条，其余单块范围为 {min(x['homopolymer_ge8'] for x in others)}–{max(x['homopolymer_ge8'] for x in others)} 条。k=19 的第12块内部交叉链接为3，其余单块为1–7。以上粗粒度特征没有显示第12块特异的明显异常；这不足以排除更细的序列依赖合成偏差。",'',
        '**后续重建的依据**','',
        '1. 优先比较保留原始长度、取消55–70 nt筛选的重建输入；当前覆盖对照已证明它能提供额外正确 k-mer 支持，但尚未对这套输入重跑 DBGPS2＋RS。',
        '2. 对第12块试验较小 k 时，同时处理候选爆炸；按块采用不同 k/c/r、混合 k 或分阶段约束搜索，比只放宽 r 更值得验证。',
        '3. 若尝试保留单次 k-mer 或扩大1050路径上限，必须同时报告错误候选、RS误纠正和最终唯一输出；覆盖改善不能直接当作恢复成功。',
        '   当前实验适配器会过滤单次键，并只枚举计数>1的种子；试验 c=0 时，必须同步调整计数和种子规则，仅更改阈值参数不足以保留完整的单次支持路径。',
        '4. 当前没有完成第12块或整个File3的新增解码验证；本次产出是实际工具诊断和与既有恢复结果的逐索引核对。','',
        f"![覆盖、k、缺失位置与搜索损失]({OUT/'diagnostics.png'})",'',
        '**数据与复核**','',
        f"- [全部输入、k、12块指标]({OUT/'all_blocks.tsv'})；[第12块逐索引诊断]({OUT/'block12_per_index.tsv'})。",
        f"- [逐位置覆盖]({OUT/'position_profiles.k19.tsv'})；[最佳参数逐块损失分解]({OUT/'best_parameter_loss_breakdown.tsv'})；[跨参数损失分解]({OUT/'union_loss_breakdown.tsv'})。",
        f"- [完整统计]({OUT/'summary.json'})；[工具来源与编译信息]({OUT/'provenance.json'})；[验证与SHA-256]({OUT/'validation.json'})。",'',
        '验证：9个k的全196596条链最小覆盖均与原DB​​GPS2过滤计数器一致（单次键在旧过滤器中记为0）；全部179组参数的完整真值路径可达性，加上原始overflow状态，逐索引精确等于实际未经修正正确候选的存在标记。原始reads计数逐位置不低于对应处理后计数。CLI输出Sm的计数与展开后的覆盖矩阵逐项一致。']
    (OUT/'REPORT.zh-CN.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
