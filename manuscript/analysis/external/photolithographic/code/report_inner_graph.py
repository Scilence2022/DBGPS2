#!/usr/bin/env python3
"""Score actual one-output-per-index decisions and label candidate upper bounds."""
import argparse
import csv
import hashlib
import json
import pathlib
import re
import sys
import zipfile
import subprocess

ROOT=pathlib.Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/inner_rs'
PAT=re.compile(r'file3\.k(\d+)\.c(\d+)\.r(\d+\.\d+)\.inner_stats\.json')
ALPHABET=str.maketrans('ACGT','0123')

def read_truth():
    return [s for s in (DATA/'sources/File3_synthesized.fasta').read_text().splitlines() if s and not s.startswith('>')]

def outputs(path):
    result={}
    for line in path.read_text().splitlines():
        i,s=line.split('\t');i=int(i)
        assert i not in result and len(s)==60
        result[i]=s
    return result

def score(output,target):
    block_size=16383;n=len(target);blocks=[dict(block=i,emitted=0,correct=0,errors=0,lane_errors=[0]*6) for i in range(n//block_size)]
    for i,s in output.items():
        assert 0<=i<n
        b=blocks[i//block_size];b['emitted']+=1
        if s==target[i]:b['correct']+=1
        else:
            b['errors']+=1
            # The 42-base payload is the little-endian stack of six 14-bit
            # field elements; randomization is XOR and preserves differences.
            delta=int(s[:42].translate(ALPHABET),4)^int(target[i][:42].translate(ALPHABET),4)
            assert delta, 'Different legal codewords with the same index must differ in payload'
            for j in range(6):b['lane_errors'][j]+=bool((delta>>(14*j))&16383)
    for b in blocks:
        b['erasures']=block_size-b['emitted']
        b['conservative_budget']=2*b['errors']+b['erasures']
        b['within_conservative_bound']=b['conservative_budget']<=5406
        b['lane_budgets']=[2*e+b['erasures'] for e in b['lane_errors']]
        b['within_rs_bound']=all(x<=5406 for x in b['lane_budgets'])
    correct=sum(b['correct'] for b in blocks);errors=len(output)-correct
    return dict(emitted=len(output),correct=correct,errors=errors,erasures=n-len(output),
                exact_recall=correct/n,precision=correct/len(output) if output else None,
                blocks_in_bound=sum(b['within_rs_bound'] for b in blocks),per_block=blocks)

def rank(row):return row['correct'],-row['errors']

def collect(target):
    rows=[];coverage=[]
    for path in sorted(OUT.glob('file3.k*.inner_stats.json')):
        m=PAT.fullmatch(path.name)
        if not m:continue
        k,c,r=int(m[1]),int(m[2]),float(m[3]);stem=path.name.removesuffix('.inner_stats.json')
        diag_path=OUT/(stem+'.truth_candidates.tsv')
        if not diag_path.exists():continue
        with diag_path.open() as handle:diag=list(csv.DictReader(handle,delimiter='\t'))
        if len(diag)!=len(target):continue
        stats=json.loads(path.read_text())
        assert stats['complete_paths']==stats['zero_error_paths']+stats['one_symbol_paths']+stats['uncorrectable_paths']+stats['outside_index_domain']
        assert sum(stats['one_symbol_position_counts'])==stats['one_symbol_paths']
        assert sum(stats['one_symbol_position_counts'][2:6])==stats['index_reassigned_paths']
        assert stats['distinct_decoded_codewords']<=stats['zero_error_paths']+stats['one_symbol_paths']
        any_correct=sum(int(x['correct_candidate_present']) for x in diag)
        zero_correct=sum(int(x['correct_zero_paths'])>0 for x in diag)
        one_only=sum(int(x['correct_one_paths'])>0 and int(x['correct_zero_paths'])==0 for x in diag)
        assert any_correct==zero_correct+one_only and any_correct<=len(target)
        coverage.append(dict(k=k,cutoff=c,ratio=r,correct_candidate_any=any_correct,
                             correct_candidate_zero=zero_correct,correct_candidate_new_from_rs=one_only,
                             exact_recall=any_correct/len(target),
                             other_distinct_codeword_candidates=stats['distinct_decoded_codewords']-any_correct,
                             per_block=[dict(block=i//16383,correct_any=sum(int(x['correct_candidate_present']) for x in diag[i:i+16383]),
                                             correct_zero=sum(int(x['correct_zero_paths'])>0 for x in diag[i:i+16383]),
                                             new_from_rs=sum(int(x['correct_one_paths'])>0 and int(x['correct_zero_paths'])==0 for x in diag[i:i+16383]))
                                        for i in range(0,len(target),16383)],stats=stats))
        for policy,emitted in stats['emitted'].items():
            p=OUT/(stem+'.'+policy+'.strands');result=score(outputs(p),target)
            assert result['emitted']==emitted and result['correct']<=any_correct
            rows.append(dict(k=k,cutoff=c,ratio=r,policy=policy,path=str(p),
                             correct_candidate_any=any_correct,correct_candidate_zero=zero_correct,
                             correct_candidate_new_from_rs=one_only,**result))
    return rows,coverage

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partial',action='store_true');ap.add_argument('--plot',action='store_true');args=ap.parse_args()
    target=read_truth();rows,coverage=collect(target)
    if not rows:print('No completed configurations');return
    best=max(rows,key=rank);block_best=max(rows,key=lambda r:(r['blocks_in_bound'],)+rank(r))
    record=dict(completed_configurations=len(coverage),evaluated_decision_policies=len(rows),
                metric_definitions=dict(candidate_coverage='Distinct reference indices whose exact 60-nt sequence occurs among decoded candidates',
                                        unique_output_correct='Exact 60-nt matches after retaining at most one sequence per decoded index',
                                        votes='Counts of graph paths, not sequencing reads'),
                best_single_configuration=best,best_block_feasibility=block_best,
                best_candidate_coverage_before_rs=max(coverage,key=lambda r:r['correct_candidate_zero']),
                candidate_coverage_upper_bound=max(coverage,key=lambda r:r['correct_candidate_any']))
    if args.partial:
        print(json.dumps(record,indent=2));return
    bounds=json.loads((OUT/'file3_candidate_path_bounds.json').read_text())
    tested={(r['k'],r['cutoff'],r['ratio']) for r in coverage}
    expected={(r['k'],r['cutoff'],r['ratio']) for r in bounds if not r['dominated_by_strict_baseline']}
    assert tested==expected,(len(tested),len(expected),expected-tested)
    assert best['correct']>=137662
    rejected=[r for r in bounds if r['dominated_by_strict_baseline']]
    assert all(r['complete_paths']<best['correct'] for r in rejected)
    # The lower bound used to prune is attained by the unchanged strict policy.
    lower=[r for r in rows if r['k']==19 and r['cutoff']==1 and r['ratio']==3.5 and r['policy']=='strict'][0]
    assert lower['correct']==137662
    record.update(total_original_configurations=len(bounds),pruned_configurations=len(rejected),
                  pruning_lower_bound=137662,largest_pruned_upper_bound=max(r['complete_paths'] for r in rejected),
                  strict_control=lower)
    baseline=outputs(DATA/'results/file3.k19.c1.r3.50.strands');winner=outputs(pathlib.Path(best['path']))
    before={i for i,s in baseline.items() if s==target[i]};after={i for i,s in winner.items() if s==target[i]}
    record['gain_vs_previous_best']=dict(gained_correct=len(after-before),lost_correct=len(before-after),net_correct=len(after)-len(before))
    strict_this=next(r for r in rows if all(r[key]==best[key] for key in ['k','cutoff','ratio']) and r['policy']=='strict')
    record['strict_same_configuration']=strict_this
    same=outputs(pathlib.Path(strict_this['path']));same_correct={i for i,s in same.items() if s==target[i]}
    diag_path=pathlib.Path(best['path'].split('.'+best['policy']+'.strands')[0]+'.truth_candidates.tsv')
    with diag_path.open() as handle:diag={int(r['id']):r for r in csv.DictReader(handle,delimiter='\t')}
    new=after-same_correct
    record['same_configuration_breakdown']=dict(strict_correct=strict_this['correct'],
        newly_correct=len(new),previous_correct_lost=len(same_correct-after),
        newly_correct_with_one_symbol_required=sum(int(diag[i]['correct_zero_paths'])==0 for i in new),
        newly_selected_existing_zero_error_candidate=sum(int(diag[i]['correct_zero_paths'])>0 for i in new),
        wrong_strict_outputs_now_correct=sum(i in same and same[i]!=target[i] for i in new),
        final_correct=best['correct'])
    record['per_policy_best']={p:max((r for r in rows if r['policy']==p),key=rank) for p in sorted({r['policy'] for r in rows})}
    eligible=[r for r in rows if r['blocks_in_bound']==12]
    if eligible:
        choice=max(eligible,key=rank);dest=OUT/'file3_recovered.zip';log=OUT/'file3_decode.log'
        with log.open('w') as handle:p=subprocess.run([str(DATA/'bin/rs_driver_fast'),'decode',choice['path'],str(dest),'12'],stdout=handle,stderr=handle,timeout=600)
        result=dict(attempted=True,exit_code=p.returncode,parameters={k:choice[k] for k in ['k','cutoff','ratio','policy']})
        if not p.returncode:
            data=dest.read_bytes();result.update(byte_identical=data==(DATA/'results/reference_file3.zip').read_bytes(),zip_error=zipfile.ZipFile(dest).testzip(),sha256=hashlib.sha256(data).hexdigest())
            assert result['byte_identical'] and result['zip_error'] is None
        record['file_decode']=result
    else:record['file_decode']=dict(attempted=False,reason='no decision output has all 12 blocks within the exact six-lane RS bound')
    (OUT/'summary.json').write_text(json.dumps(record,indent=2))
    (OUT/'all_metrics.json').write_text(json.dumps(rows,indent=2))
    (OUT/'candidate_coverage.json').write_text(json.dumps(coverage,indent=2))
    columns=['k','cutoff','ratio','policy','emitted','correct','errors','erasures','exact_recall','precision','blocks_in_bound','correct_candidate_any','correct_candidate_zero','correct_candidate_new_from_rs','path']
    with (OUT/'comparison.csv').open('w') as handle:
        w=csv.DictWriter(handle,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)
    if args.plot:plot(record)
    write_report(record)
    print(json.dumps(record,indent=2))

def write_report(record):
    b=record['best_single_configuration'];s=record['strict_same_configuration'];old=record['strict_control'];p=record['best_block_feasibility'];upper=record['candidate_coverage_upper_bound']
    text=['**File 3：DBGPS2 图候选＋内层 RS 单符号纠错**','',
          f"在既有参数网格内，按论文“至少出现一次”口径的最高正确片段覆盖为 **{upper['correct_candidate_any']:,}/196596（{100*upper['exact_recall']:.6f}%）**；要求每个索引最终只留一条时，最高正确输出为 **{b['correct']:,}/196596（{100*b['exact_recall']:.6f}%）**。两项最优参数分别列于下文，均按索引去重。",'',
          '本报告区分“内码校验通过（两个校验子均为零）”与“修正一个 GF(64) 符号”。此前 137662 条正确片段属于前者。RS(20,18) 的两个校验符号可纠正一个 6-bit 符号错误；该符号覆盖 3 nt，可以包含多处碱基替换，不代表能任意纠正一次插入或缺失。','',
          '解码计算为 `S0=Σc_j`、`S1=Σα^j·c_j`，运算在 GF(64) 中进行，α=2，本原多项式为 91。两项均为 0 时不修改；否则寻找 `S1=S0·α^j` 的 j∈[0,19]，并执行 `c_j←c_j XOR S0`。不存在该位置时拒绝。多符号错误可能被误纠正到其他合法码字，实际正确数必须用真值另行核对。','',
          f"保留原先 179 组 k/c/r 参数、1050 路径上限和完整保留 reads。{record['completed_configurations']} 组实际重跑，其余 {record['pruned_configurations']} 组的完整候选路径数已少于可达到的严格基线 137662，故不可能刷新正确数；被排除组合的最大候选数上界为 {record['largest_pruned_upper_bound']}。每条完整路径在半径 1 内至多对应一个 RS 码字。溢出搜索不提供候选，沿用原规则。",'',
          '该剪枝证明的是最高完整正确片段数；它不用于判断被排除配置的误片段分布或完整文件可恢复性。外码达标和文件解码只针对实际重跑的配置报告。','',
          '对每个完整候选计算两个校验子，执行确定的有界距离解码；索引若被纠正，按纠正后的索引重新归组。原始候选与解码码字不借助参考选择。参考只在所有决策输出关闭后加载，用于评估。所有输出均与正确的索引及完整 60 nt 真值比较，并对索引去重。','',
          '比较的决策规则包括：优先保留唯一的无需修正码字；没有无需修正的候选时接受唯一纠正结果，或按图路径支持票数投票（最低票数 1/2/3/5/10/20，并列缺失）；以及把无需修正与纠正候选共同投票，前者的票权为 1/4/16。这里的票数是图路径数，不是独立测序 reads 数。多错误候选可能被误纠正，故必须同时报告错误输出。','',
          '| 方法 | k / c / r | 输出链数 | 正确链数 | 错误链数 | 缺失 | 正确恢复率 |',
          '|---|---|---:|---:|---:|---:|---:|']
    for label,r in [('此前最高（只校验）',old),('新最佳同参数只校验',s),(b['policy'],b)]:
        text.append(f"| {label} | {r['k']} / {r['cutoff']} / {r['ratio']:.2f} | {r['emitted']:,} | {r['correct']:,} | {r['errors']:,} | {r['erasures']:,} | {100*r['exact_recall']:.4f}% |")
    gain=record['gain_vs_previous_best'];text+=['',f"相对此前最高，新增正确链 {gain['gained_correct']:,} 条，丢失此前正确链 {gain['lost_correct']:,} 条，净增 {gain['net_correct']:,} 条。最高数是本参数与决策网格中的实际单配置输出，不是逐索引用真值挑选候选。",'',
             f"在最高正确数的同一组图参数下，纠错后正确码字至少存在于候选池的索引有 {b['correct_candidate_any']:,} 个；最终正确输出 {b['correct']:,} 个，相差 {b['correct_candidate_any']-b['correct']:,} 个，属于已有正确候选但未被该决策选出的情况。",'',
             '**各决策规则的已测最高正确数**','',
             '| 规则 | k / c / r | 正确 | 错误 | 恢复率 | 达标外码块 |',
             '|---|---|---:|---:|---:|---:|']
    d=record['same_configuration_breakdown']
    text[-4:-4]=['',f"同参数分解：`{d['strict_correct']} + {d['newly_correct']} - {d['previous_correct_lost']} = {d['final_correct']}` 条正确链。新增正确链中，{d['newly_correct_with_one_symbol_required']} 条必须经过单符号纠错才出现在候选中，{d['newly_selected_existing_zero_error_candidate']} 条来自重新选择原先已有的零错误候选。",'',
                  f"数量核算：`输出 {b['emitted']} = 正确 {b['correct']} + 错误 {b['errors']}`；`总数 196596 = 输出 {b['emitted']} + 缺失 {b['erasures']}`。正确恢复率为 `{b['correct']}/196596 = {100*b['exact_recall']:.6f}%`，输出精确率为 `{b['correct']}/{b['emitted']} = {100*b['precision']:.6f}%`。",'']
    for policy,r in record['per_policy_best'].items():text.append(f"| {policy} | {r['k']} / {r['cutoff']} / {r['ratio']:.2f} | {r['correct']:,} | {r['errors']:,} | {100*r['exact_recall']:.4f}% | {r['blocks_in_bound']}/12 |")
    if (OUT/'comparison.png').exists():text+=['',f"![正确数、误纠正与逐块纠错预算]({OUT/'comparison.png'})"]
    before=record['best_candidate_coverage_before_rs']
    text+=['', '**按论文“至少出现一次”口径的最高正确片段覆盖**','',
           f"不做符号修正时，原始完整候选中已精确匹配真值的最高覆盖为 {before['correct_candidate_zero']:,} 条，对应 k={before['k']}、c={before['cutoff']}、r={before['ratio']:.2f}。启用单符号纠错后，整个参数网格的最高候选覆盖净增 {upper['correct_candidate_any']-before['correct_candidate_zero']:,} 条。该比较与每索引唯一输出的 137662 条基线具有不同的候选保留口径。",'',
           f"最高为 k={upper['k']}、c={upper['cutoff']}、r={upper['ratio']:.2f}，纠错后候选池至少包含正确序列的索引有 **{upper['correct_candidate_any']:,}/196596（{100*upper['exact_recall']:.6f}%）**。同参数分解为 `原先已有正确完整路径 {upper['correct_candidate_zero']} + 仅在单符号纠错后出现 {upper['correct_candidate_new_from_rs']} = {upper['correct_candidate_any']}`。",'',
           f"这是实际生成候选的正确覆盖数，与论文“原始片段至少出现一次”口径对应；它同时是该候选池每索引唯一决策的正确数上限。该候选池另含 {upper['other_distinct_codeword_candidates']:,} 个不同的非真值码字，不能依据未知真值直接删除它们。因此候选覆盖不能代替最终唯一输出的正确数或完整文件恢复。",'',
           '这里对齐的是“至少出现一次”的计数定义；本表的处理阶段是内码纠错之后。与论文聚类阶段的恢复率比较时，还需要统一处理阶段。','',
           '**最高候选覆盖配置的逐阶段数量**','',
           '| 阶段 | 数量 |', '|---|---:|']
    st=upper['stats']
    for label,value in [('完整图路径（未去重）',st['complete_paths']),('两个校验子均为零且索引有效的路径',st['zero_error_paths']),('修正一个符号且索引有效的路径',st['one_symbol_paths']),('无法按单符号解码的路径',st['uncorrectable_paths']),('解码后索引超出有效范围的路径',st['outside_index_domain']),('去重后的全部合法码字',st['distinct_decoded_codewords']),('至少出现一次的正确码字',upper['correct_candidate_any']),('其他非真值码字',upper['other_distinct_codeword_candidates'])]:
        text.append(f'| {label} | {value:,} |')
    text+=['',f"路径数核对：`{st['complete_paths']} = {st['zero_error_paths']} + {st['one_symbol_paths']} + {st['uncorrectable_paths']} + {st['outside_index_domain']}`。有效解码路径合计 {st['zero_error_paths']+st['one_symbol_paths']:,} 条，按索引和完整 DNA 去重后为 {st['distinct_decoded_codewords']:,} 个码字。单符号修正路径中，数据区占 {sum(st['one_symbol_position_counts'][6:]):,} 条、索引区占 {sum(st['one_symbol_position_counts'][2:6]):,} 条、校验区占 {sum(st['one_symbol_position_counts'][:2]):,} 条；这些计数包含误纠正与重复路径。",'',
           '**最高唯一输出正确数配置的逐块结果**','',
           '| 块 | 正确 | 错误片段 | 缺失 | 六个外码列中最大 2E+S | 达标 |',
           '|---:|---:|---:|---:|---:|---|']
    for x in b['per_block']:text.append(f"| {x['block']+1} | {x['correct']:,} | {x['errors']:,} | {x['erasures']:,} | {max(x['lane_budgets']):,} | {'是' if x['within_rs_bound'] else '否'} |")
    text+=['',f"实际重跑组中，外码达标数最多的配置为 {p['k']}/{p['cutoff']}/{p['ratio']:.2f}/{p['policy']}，{p['blocks_in_bound']}/12 块。每块含六个独立 GF(2^14) 码字，预算按各列真实错误数分别计算，阈值均为 5406；完整片段错误数仅给出更保守的共同上界。",'',
           '图搜索的严格输出逐配置与原结果做字节比较；单符号解码在 100 条公开码字的 20 个位置、全部 63 种非零错误上通过 126000 次精确恢复检查。所有最佳参数均按同一数据真值事后选出，属于探索结果，没有独立留出集。','',
           f"完整指标：[JSON]({OUT/'summary.json'})；所有参数与决策：[CSV]({OUT/'comparison.csv'})；候选上界证明：[原参数候选数]({OUT/'file3_candidate_path_bounds.json'})；原始设计与测序数据保留在 `{DATA}`。未修改生产 DBGPS2 或原严格适配程序，未提交或推送。"]
    if record['file_decode'].get('byte_identical'):text+=['',f"完整 File 3 已通过 ZIP CRC 和参考文件逐字节校验：[ZIP]({OUT/'file3_recovered.zip'})。"]
    else:text+=['','本扫描未取得通过验证的完整 File 3；片段最高正确数与完整文件恢复分别报告。']
    (OUT/'REPORT.zh-CN.md').write_text('\n'.join(text)+'\n')

def plot(record):
    sys.path.insert(0,str(DATA/'plotlibs'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    old=record['strict_control'];best=record['best_single_configuration'];block=record['best_block_feasibility']
    cases=[old,record['per_policy_best']['corrected_unique'],record['per_policy_best']['corrected_vote2'],best]
    labels=['Previous strict best','Unique RS repair','Priority vote >=2','Maximum correct']
    fig,axes=plt.subplots(1,2,figsize=(14,4.8),layout='constrained')
    for i,(label,r) in enumerate(zip(labels,cases)):
        offset=0
        for field,color in [('correct','#317c65'),('errors','#d47d38'),('erasures','#c9ced4')]:
            axes[0].barh(i,r[field],left=offset,color=color,label=field.title() if i==0 else None)
            offset+=r[field]
        axes[0].text(r['correct']/2,i,f"{r['correct']:,}",ha='center',va='center',color='white',fontsize=9)
        axes[0].text(198500,i,f"E={r['errors']:,}",va='center',fontsize=9)
    axes[0].set(yticks=range(4),yticklabels=labels,xlabel='Indexed 60-nt strands',xlim=(0,235000),title='Correct, incorrect, and missing outputs')
    axes[0].invert_yaxis();axes[0].legend(frameon=False,loc='upper center',bbox_to_anchor=(.5,-.12),ncol=3)
    for r,label,color in [(old,'Previous strict best','#6b7785'),(best,'Maximum correct','#d47d38'),(block,'Most blocks within bound','#317c65')]:
        axes[1].plot(range(1,13),[max(x['lane_budgets']) for x in r['per_block']],'-o',label=label,color=color)
    axes[1].axhline(5406,color='#b33e40',linestyle='--',label='RS limit = 5406')
    axes[1].set(xlabel='Outer block',ylabel='Worst of six lanes: 2E + S',xticks=range(1,13),title='More correct strands can still mean more errors')
    axes[1].grid(alpha=.2);axes[1].legend(frameon=False,fontsize=8)
    fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.svg');plt.close(fig)

if __name__=='__main__':main()
