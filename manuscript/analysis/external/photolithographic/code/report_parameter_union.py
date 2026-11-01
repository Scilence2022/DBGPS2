#!/usr/bin/env python3
"""Audit reference-labelled, per-decoded-index coverage across graph settings."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
SOURCE=DATA/'results/inner_rs'
N=196596
BLOCK=16383
popcount=getattr(int,'bit_count',lambda value:bin(value).count('1'))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--partial',action='store_true')
    args=ap.parse_args()
    out=SOURCE/('parameter_union_partial' if args.partial else 'parameter_union')
    out.mkdir(exist_ok=True)
    plan=sorted(json.loads((SOURCE/'file3_candidate_path_bounds.json').read_text()),key=lambda r:(r['k'],r['cutoff'],r['ratio']))
    assert len({(r['k'],r['cutoff'],r['ratio']) for r in plan})==len(plan)==179
    truth=[s for s in (DATA/'sources/File3_synthesized.fasta').read_text().splitlines() if s and not s.startswith('>')]
    assert len(truth)==N and all(len(s)==60 for s in truth)
    candidate=[0]*N;zero=[0]*N;policies={};combined=[0]*N
    configurations=[];candidate_bits={};original60=0;new119=0;validated=0
    with gzip.open(out/'successful_parameter_pairs.tsv.gz','wt') as pairs:
        pairs.write('decoded_index\tparameter_id\tcorrect_zero_paths\tcorrect_one_paths\n')
        for p,row in enumerate(plan):
            name=f"file3.k{row['k']}.c{row['cutoff']}.r{row['ratio']:.2f}"
            stats_path=SOURCE/(name+'.inner_stats.json');diag_path=SOURCE/(name+'.truth_candidates.tsv')
            if not stats_path.exists() or not diag_path.exists():
                if args.partial:continue
                raise AssertionError('Missing configuration: '+name)
            # Snapshot a finished diagnostic; partial reports exclude in-flight files.
            with diag_path.open('rb') as handle:
                handle.seek(max(0,diag_path.stat().st_size-256));lines=handle.read().splitlines()
            last=lines[-1].split(b'\t') if lines else []
            if len(last)!=4 or last[0]!=b'196595':
                if args.partial:continue
                raise AssertionError('Incomplete diagnostic: '+name)
            stats=json.loads(stats_path.read_text());bit=1<<p
            zcount=0;ccount=0;bits=bytearray((N+7)//8)
            with diag_path.open() as handle:
                assert next(handle).strip()=='id\tcorrect_zero_paths\tcorrect_one_paths\tcorrect_candidate_present'
                seen=0
                for expected,line in enumerate(handle):
                    i,z,o,c=map(int,line.split());assert i==expected and i<N and c==bool(z or o)
                    seen+=1
                    if z:zero[i]|=bit;zcount+=1
                    if c:
                        candidate[i]|=bit;ccount+=1;bits[i>>3]|=1<<(i&7)
                        pairs.write(f'{i}\t{p}\t{z}\t{o}\n')
            assert seen==N
            assert stats['complete_paths']==sum(stats[k] for k in ['zero_error_paths','one_symbol_paths','uncorrectable_paths','outside_index_domain'])
            cb=int.from_bytes(bits,'little');candidate_bits[p]=cb
            if row['dominated_by_strict_baseline']:new119|=cb
            else:original60|=cb
            policy_counts={}
            for policy,emitted in stats['emitted'].items():
                success=policies.setdefault(policy,[0]*N);correct=0;seen_ids=set()
                path=SOURCE/(name+'.'+policy+'.strands')
                with path.open() as handle:
                    for line in handle:
                        field,s=line.rstrip('\n').split('\t');i=int(field)
                        assert 0<=i<N and i not in seen_ids and len(s)==60
                        seen_ids.add(i)
                        if s==truth[i]:
                            assert candidate[i]&bit
                            if policy=='strict':assert zero[i]&bit
                            success[i]|=bit;combined[i]|=bit;correct+=1
                assert len(seen_ids)==emitted
                policy_counts[policy]=dict(correct=correct,errors=emitted-correct,emitted=emitted)
            control=(SOURCE/(name+'.strands')).read_bytes()
            assert control==(SOURCE/(name+'.strict.strands')).read_bytes()==(DATA/'results'/(name+'.strands')).read_bytes()
            validated+=1
            configurations.append(dict(parameter_id=p,k=row['k'],cutoff=row['cutoff'],ratio=row['ratio'],
                previously_pruned=row['dominated_by_strict_baseline'],candidate_correct=ccount,
                candidate_per_block=[popcount((cb>>(b*BLOCK))&((1<<BLOCK)-1)) for b in range(12)],
                zero_candidate_correct=zcount,new_from_rs=ccount-zcount,policies=policy_counts,
                strict_sha256=hashlib.sha256(control).hexdigest()))
            print('AUDITED',name,'candidate',ccount,'configurations',validated,flush=True)
    if not args.partial:assert len(configurations)==len(plan)==179
    assert policies
    def total(values):return sum(bool(x) for x in values)
    def blocks(values):return [sum(bool(x) for x in values[a:a+BLOCK]) for a in range(0,N,BLOCK)]
    candidate_count=total(candidate);zero_count=total(zero)
    policy_summary={policy:dict(correct_union=total(values),recall=total(values)/N,per_block=blocks(values),
                               best_single=max((r['policies'][policy]['correct'] for r in configurations)))
                    for policy,values in sorted(policies.items())}
    best_policy=max(policy_summary,key=lambda p:policy_summary[p]['correct_union'])
    best_single=max(configurations,key=lambda r:r['candidate_correct'])
    all_bits=0
    for bits in candidate_bits.values():all_bits|=bits
    assert popcount(all_bits)==candidate_count
    covered=0;remaining=dict(candidate_bits);greedy=[]
    while remaining:
        p=max(remaining,key=lambda p:(popcount(remaining[p]&~covered),-p))
        gained=popcount(remaining[p]&~covered)
        if not gained:break
        covered|=remaining.pop(p)
        greedy.append(dict(parameter_id=p,gained=gained,cumulative=popcount(covered)))
    assert covered==all_bits
    per_k=[]
    for k in sorted({r['k'] for r in configurations}):
        mask=sum(1<<r['parameter_id'] for r in configurations if r['k']==k)
        per_k.append(dict(k=k,correct_union=sum(bool(x&mask) for x in candidate),
                          exclusive_to_this_k=sum(bool(x&mask) and not(x&~mask) for x in candidate)))
    completed_mask=sum(1<<r['parameter_id'] for r in configurations)
    exclusive=Counter(x.bit_length()-1 for x in candidate if x and popcount(x)==1)
    for r in configurations:r['exclusive_correct_indices']=exclusive[r['parameter_id']]
    summary=dict(total_indices=N,total_planned_parameters=len(plan),completed_parameters=len(configurations),
        candidate_union_correct=candidate_count,candidate_union_recall=candidate_count/N,
        zero_candidate_union_correct=zero_count,requires_rs_across_all_parameters=sum(bool(c) and not z for c,z in zip(candidate,zero)),
        missing_after_candidate_union=N-candidate_count,
        best_single_candidate={k:best_single[k] for k in ['parameter_id','k','cutoff','ratio','candidate_correct']},
        gain_over_best_single_candidate=candidate_count-best_single['candidate_correct'],
        previous_60_union=popcount(original60),new_119_union=popcount(new119),
        gained_only_from_previously_pruned=popcount(new119&~original60),
        per_policy=policy_summary,best_fixed_policy=best_policy,all_policy_output_union_correct=total(combined),
        per_block=[dict(block=b+1,correct=c,missing=BLOCK-c,oracle_erasure_bound_met=BLOCK-c<=5406,
                       best_single_configuration_correct=best_single['candidate_per_block'][b],
                       gain_over_best_single_configuration=c-best_single['candidate_per_block'][b])
                   for b,c in enumerate(blocks(candidate))],per_k=per_k,
        success_parameter_count_histogram=dict(sorted(Counter(popcount(x) for x in candidate).items())),
        greedy_cover=greedy,strict_byte_controls=validated,
        definitions=dict(index='Corrected decoded index, zero-based, rather than the original graph seed index',
            candidate_success='At least one complete non-overflow path becomes the exact full 60-nt reference strand after radius-one inner RS decoding',
            output_success='At least one setting emits the exact 60-nt reference at this decoded index under the named fixed decision policy',
            matrix='0: no exact candidate; 1: exact candidate exists without correction; 2: exact candidate exists only after one-symbol correction',
            oracle='Reference-labelled union coverage; a blind algorithm cannot choose the correct candidate or parameter using these labels',
            greedy='Deterministic greedy coverage ordering, not a proven minimum parameter set'))
    with (out/'per_index.tsv').open('w') as handle:
        names=['decoded_index','block_1based','offset_0based','candidate_success','zero_candidate_success','rs_required_globally','successful_candidate_parameters','successful_strict_parameters','successful_best_fixed_policy_parameters','successful_any_policy_parameters','first_candidate_parameter','first_best_fixed_policy_parameter']
        handle.write('\t'.join(names)+'\n')
        for i,(c,z,u) in enumerate(zip(candidate,zero,combined)):
            fixed=policies[best_policy][i]
            first=lambda x:(x&-x).bit_length()-1 if x else -1
            fields=[i,i//BLOCK+1,i%BLOCK,int(bool(c)),int(bool(z)),int(bool(c) and not z),popcount(c),popcount(policies['strict'][i]),popcount(fixed),popcount(u),first(c),first(fixed)]
            handle.write('\t'.join(map(str,fields))+'\n')
    cols=[r['parameter_id'] for r in configurations]
    with gzip.open(out/'candidate_success_matrix.tsv.gz','wt') as handle:
        handle.write('decoded_index\t'+'\t'.join(f'p{p}' for p in cols)+'\n')
        for i,(c,z) in enumerate(zip(candidate,zero)):
            handle.write(str(i)+'\t'+'\t'.join('1' if z&(1<<p) else '2' if c&(1<<p) else '0' for p in cols)+'\n')
    with gzip.open(out/'policy_success_masks.tsv.gz','wt') as handle:
        names=sorted(policies);handle.write('decoded_index\t'+'\t'.join(names)+'\n')
        for i in range(N):handle.write(str(i)+'\t'+'\t'.join(hex(policies[p][i]) for p in names)+'\n')
    with (out/'parameters.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=['parameter_id','k','cutoff','ratio','candidate_correct','zero_candidate_correct','new_from_rs','exclusive_correct_indices','previously_pruned'],delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(configurations)
    (out/'never_recovered_indices.txt').write_text(''.join(f'{i}\n' for i,c in enumerate(candidate) if not c))
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    (out/'configurations.json').write_text(json.dumps(configurations,indent=2))
    lines=['**File 3 跨参数逐索引恢复并集**','',
        f"核对 {len(configurations)}/{len(plan)} 组 k/c/r 参数，按内码纠正后的索引去重：候选池正确恢复并集为 **{candidate_count:,}/{N:,}（{100*candidate_count/N:.6f}%）**，仍缺失 {N-candidate_count:,} 条。",'',
        f"计数定义为 `Σ_i 1[存在参数 θ，使候选经内码 RS 解码后等于第 i 条完整 60 nt 真值]`。同一索引出现多少次均只计一次。索引也可能被 RS 修正，因此以修正后的 decoded index 归组，而非原始图搜索 seed。",'',
        f"最佳单参数候选覆盖为 {best_single['candidate_correct']:,} 条（k={best_single['k']}、c={best_single['cutoff']}、r={best_single['ratio']:.2f}）；跨参数增加 **{candidate_count-best_single['candidate_correct']:,}** 条。无需符号修正的候选并集为 {zero_count:,} 条；另有 **{candidate_count-zero_count:,}** 条在所有参数中均必须依赖单符号纠正。",'',
        f"原先 60 组并集 {popcount(original60):,} 条；此前只为寻找单参数最大值而排除的 119 组贡献了 {popcount(new119&~original60):,} 条独有正确片段。本次完整并集不采用该数量上界剪枝。",'',
        '这是利用真值核对的覆盖率。若同一索引存在多个合法码字，或不同参数输出互相冲突，覆盖标签不能直接指示实际应选哪条链。没有把真值筛选得到的并集作为盲解码输出或完整文件恢复。','',
        '**每个参数最终只输出一条链时，正确输出的跨参数并集**','',
        '| 固定决策规则 | 最佳单参数正确数 | 跨参数正确并集 | 恢复率 |',
        '|---|---:|---:|---:|']
    for p,r in policy_summary.items():lines.append(f"| {p} | {r['best_single']:,} | {r['correct_union']:,} | {100*r['recall']:.4f}% |")
    lines+=['',f"进一步允许切换决策规则时，全部规则最终输出的正确片段并集为 {total(combined):,} 条。该数同样是逐索引用真值事后核对，不能将这些片段假定为可无误识别。",'',
        '**逐块候选覆盖**','', '| 块 | 最佳单参数中的正确数 | 跨参数正确并集 | 新增 | 尚缺 | 理想筛选后的外码擦除条件 |','|---:|---:|---:|---:|---:|---|']
    for b in summary['per_block']:lines.append(f"| {b['block']} | {b['best_single_configuration_correct']:,} | {b['correct']:,} | {b['gain_over_best_single_configuration']:,} | {b['missing']:,} | {'满足' if b['oracle_erasure_bound_met'] else '不满足'} |")
    lines+=['','每块 16383 条链，外码 RS(16383,10977) 可处理 5406 个擦除。此表假定能识别并排除所有错误候选，故只说明理想真值筛选条件下的覆盖潜力。真实合并还需确定候选选择规则、计算各外码列的错误数并实际解码。','',
        '**按 k 的覆盖与独有贡献**','', '| k | 该 k 下参数并集 | 只有该 k 能恢复 |','|---:|---:|---:|']
    for r in per_k:lines.append(f"| {r['k']} | {r['correct_union']:,} | {r['exclusive_to_this_k']:,} |")
    lookup={r['parameter_id']:r for r in configurations}
    lines+=['','**参数间的互补性：贪心覆盖顺序前十组**','',
        '| 顺序 | k / c / r | 新增正确索引 | 累计正确索引 |', '|---:|---|---:|---:|']
    for j,step in enumerate(greedy[:10]):
        r=lookup[step['parameter_id']]
        lines.append(f"| {j+1} | {r['k']} / {r['cutoff']} / {r['ratio']:.2f} | {step['gained']:,} | {step['cumulative']:,} |")
    lines+=['',f"每一步选择对当前并集新增正确索引最多的参数，最终 {len(greedy)} 组覆盖本扫描的全部正确并集。这使用真值标签，是对参数互补性的描述，不证明该参数子集数量最小，也不提供盲选择规则。"]
    lines+=['','**可逐索引复核的数据**','',
        f"- [逐索引汇总]({out/'per_index.tsv'})：每个索引一行，包含成功次数、是否必须依赖 RS、首个成功参数。",
        f"- [候选成功矩阵（gzip）]({out/'candidate_success_matrix.tsv.gz'})：每行一个索引，每列一组参数；0=无正确候选，1=有无需修正的正确候选，2=仅纠错后有正确候选。",
        f"- [全部成功的索引与参数配对（gzip）]({out/'successful_parameter_pairs.tsv.gz'})：还保留正确零错误路径数和正确单符号纠错路径数。",
        f"- [参数表]({out/'parameters.tsv'})：矩阵 p 列、首个成功参数及位掩码都使用该表的 parameter_id。",
        f"- [各输出规则成功位掩码（gzip）]({out/'policy_success_masks.tsv.gz'})：第 p 位为 1 表示该规则在参数 p 下最终输出正确。",
        f"- [始终未恢复的索引]({out/'never_recovered_indices.txt'})；[完整统计]({out/'summary.json'})。",'',
        f"校验：{validated} 组严格输出与原始结果及新收集器严格分支逐字节一致；每组诊断恰有 {N} 个连续索引；每个正确最终输出均属于同参数的正确候选集合。沿用原图适配器、1050 路径上限、完整 reads 和 RS 单符号解码方法。",
        '',f"[来源与文件 SHA-256]({out/'validation.json'})；[独立交叉核对结果]({out/'independent_audit.json'})（运行 `validate_parameter_union.py` 生成）。"]
    (out/'REPORT.zh-CN.md').write_text('\n'.join(lines)+'\n')
    manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.iterdir()) if p.is_file() and p.name!='validation.json'}
    source_files=[Path(__file__),ROOT/'DBGPS2.c',ROOT/'benchmarks/antkowiak2020/graph_inner_rs.cpp',ROOT/'benchmarks/antkowiak2020/graph_adapter.c',DATA/'bin/graph_inner_rs',DATA/'sources/File3_synthesized.fasta',SOURCE/'file3_candidate_path_bounds.json']
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    (out/'validation.json').write_text(json.dumps(dict(completed_mask=hex(completed_mask),strict_controls=validated,files=manifest,source_sha256=source_hashes),indent=2))
    print(json.dumps({k:summary[k] for k in ['completed_parameters','candidate_union_correct','zero_candidate_union_correct','requires_rs_across_all_parameters','gained_only_from_previously_pruned','best_fixed_policy','all_policy_output_union_correct']},indent=2),flush=True)

if __name__=='__main__':main()
