#!/usr/bin/env python3
"""Score real per-block decisions, oracle coverage, and matched mixed controls."""
import argparse
import csv
import hashlib
import itertools
import json
import re
import subprocess
import zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2];DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/block_decode';B=16383;N=196596
ALPHABET=str.maketrans('ACGT','0123')
COMMON={(15,1,1.8),(15,1,2.55),(17,1,2.55),(17,1,3.5),(18,1,2.55),(18,1,3.5),(19,1,2.55),(19,1,3.5),(20,1,3.5),(20,1,5.0)}
PAT=re.compile(r'file3\.k(\d+)\.c(\d+)\.r([\d.]+)\.inner_stats\.json')

def bits(flags):return int.from_bytes(np.packbits(np.asarray(flags,dtype=np.uint8),bitorder='little').tobytes(),'little')
def nbits(value):return value.bit_count()

def score_outputs(path,truth,blocks):
    rows={b:dict(block=b,emitted=0,correct=0,errors=0,lane_errors=[0]*6) for b in blocks}
    flags={b:bytearray(B) for b in blocks};seen=set();outside=0
    with path.open() as handle:
        for line in handle:
            field,s=line.rstrip().split('\t');i=int(field)
            assert 0<=i<N and i not in seen and len(s)==60;seen.add(i)
            b=i//B+1
            if b not in rows:outside+=1;continue
            row=rows[b];row['emitted']+=1
            if s==truth[i]:row['correct']+=1;flags[b][i%B]=1
            else:
                row['errors']+=1
                delta=int(s[:42].translate(ALPHABET),4)^int(truth[i][:42].translate(ALPHABET),4);assert delta
                for j in range(6):row['lane_errors'][j]+=bool((delta>>(14*j))&16383)
    for b,r in rows.items():
        r['erasures']=B-r['emitted'];r['lane_budgets']=[2*e+r['erasures'] for e in r['lane_errors']]
        r['max_lane_budget']=max(r['lane_budgets']);r['within_rs_bound']=r['max_lane_budget']<=5406
        r['precision']=r['correct']/r['emitted'] if r['emitted'] else None;r['recall']=r['correct']/B
        r['outside_target_block_outputs']=outside
    return rows,{b:bits(np.frombuffer(v,dtype=np.uint8)) for b,v in flags.items()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partial',action='store_true');args=ap.parse_args()
    truth=[s for s in (DATA/'sources/File3_synthesized.fasta').read_text().splitlines() if s and not s.startswith('>')];assert len(truth)==N
    configs=[];decisions=[];groups={};done=[]
    def aggregate(mode,b,subset,param,cb,zb,policy_rows,policy_bits):
        key=(mode,b,subset)
        if key not in groups:groups[key]=dict(mode=mode,block=b,subset=subset,parameters=set(),candidate=0,zero=0,policies={},best_candidate=None,best_unique=None,best_rs=None)
        g=groups[key];g['parameters'].add(tuple(param.values()));g['candidate']|=cb;g['zero']|=zb
        cr=dict(**param,correct=nbits(cb),zero_correct=nbits(zb))
        if g['best_candidate'] is None or cr['correct']>g['best_candidate']['correct']:g['best_candidate']=cr
        for policy,row in policy_rows.items():
            assert not(policy_bits[policy]&~cb)
            g['policies'][policy]=g['policies'].get(policy,0)|policy_bits[policy]
            if g['best_unique'] is None or (row['correct'],-row['errors'])>(g['best_unique']['correct'],-g['best_unique']['errors']):g['best_unique']=row
            if g['best_rs'] is None or (row['max_lane_budget'],-row['correct'],row['errors'])<(g['best_rs']['max_lane_budget'],-g['best_rs']['correct'],g['best_rs']['errors']):g['best_rs']=row
    sources=[]
    for path in sorted((OUT/'decodes').glob('*/*/file3.k*.inner_stats.json')):
        mode=path.parent.parent.name;block=int(path.parent.name[5:])
        if mode=='length60_control':continue
        sources.append((path,mode,block))
    for k,c,r in sorted(COMMON):
        sources.append((DATA/f'results/inner_rs/file3.k{k}.c{c}.r{r:.2f}.inner_stats.json','historical_processed',0))
    for path,mode,block in sources:
        m=PAT.fullmatch(path.name);assert m
        k,c,r=int(m[1]),int(m[2]),float(m[3]);stem=str(path)[:-len('.inner_stats.json')]
        diag=Path(stem+'.truth_candidates.tsv')
        if not diag.exists():continue
        with diag.open('rb') as handle:
            handle.seek(max(0,diag.stat().st_size-256));last=handle.read().splitlines()
        if not last or not last[-1].startswith(b'196595\t'):continue
        stats=json.loads(path.read_text());blocks=[block] if block else list(range(1,13));cb={};zb={}
        with diag.open() as handle:
            next(handle)
            if block:rows=itertools.islice(handle,(block-1)*B,block*B)
            else:rows=handle
            flags={b:np.zeros(B,dtype=np.uint8) for b in blocks};zero={b:np.zeros(B,dtype=np.uint8) for b in blocks}
            for line in rows:
                i,z,o,a=map(int,line.split());b=i//B+1;assert a==bool(z or o)
                flags[b][i%B]=a;zero[b][i%B]=bool(z)
            cb={b:bits(v) for b,v in flags.items()};zb={b:bits(v) for b,v in zero.items()}
        policy_rows={b:{} for b in blocks};policy_bits={b:{} for b in blocks}
        for policy,emitted in stats['emitted'].items():
            output=Path(stem+'.'+policy+'.strands');scores,correct_bits=score_outputs(output,truth,blocks)
            assert sum(s['emitted'] for s in scores.values())+(next(iter(scores.values()))['outside_target_block_outputs'] if block else 0)==emitted
            for b,row in scores.items():
                result=dict(mode=mode,k=k,cutoff=c,ratio=r,policy=policy,path=str(output),**row)
                policy_rows[b][policy]=result;policy_bits[b][policy]=correct_bits[b];decisions.append(result)
        param=dict(k=k,cutoff=c,ratio=r)
        for b in blocks:
            configs.append(dict(mode=mode,block=b,**param,candidate_correct=nbits(cb[b]),zero_candidate_correct=nbits(zb[b])))
            aggregate(mode,b,'all',param,cb[b],zb[b],policy_rows[b],policy_bits[b])
            if (k,c,r) in COMMON:aggregate(mode,b,'common10',param,cb[b],zb[b],policy_rows[b],policy_bits[b])
            if mode=='broad' and b==12:aggregate(mode,b,'k'+str(k),param,cb[b],zb[b],policy_rows[b],policy_bits[b])
        done.append(dict(mode=mode,block=block,**param));print('SCORED',mode,block,k,c,r,flush=True)
    new=[r for r in done if r['mode'] not in ['historical_processed','pooled_raw']]
    ancillary=[r for r in done if r['mode']=='pooled_raw']
    if not args.partial:
        assert len(new)==319,len(new)
        assert len(groups[('broad',12,'all')]['parameters'])==179
        plan=json.loads((DATA/'results/inner_rs/file3_candidate_path_bounds.json').read_text())
        assert groups[('broad',12,'all')]['parameters']=={(r['k'],r['cutoff'],r['ratio']) for r in plan}
        for b in range(1,13):assert len(groups[('broad',b,'common10')]['parameters'])==10
        for mode in ['strict','inclusive','assigned_pooled']:assert len(groups[(mode,12,'common10')]['parameters'])==10
    summary_groups=[];raw_groups={}
    for key,g in sorted(groups.items()):
        raw_groups[key]=g
        p={name:nbits(value) for name,value in g['policies'].items()};all_policy=0
        for value in g['policies'].values():all_policy|=value
        summary_groups.append(dict(mode=g['mode'],block=g['block'],subset=g['subset'],parameters=len(g['parameters']),
            candidate_union_correct=nbits(g['candidate']),candidate_union_recall=nbits(g['candidate'])/B,
            zero_candidate_union_correct=nbits(g['zero']),rs_required_globally=nbits(g['candidate']&~g['zero']),
            policy_union_correct=p,any_policy_output_union_correct=nbits(all_policy),
            best_single_candidate=g['best_candidate'],best_single_unique=g['best_unique'],best_outer_budget=g['best_rs']))
    old_flags=np.loadtxt(DATA/'results/inner_rs/parameter_union/per_index.tsv',skiprows=1,usecols=3,dtype=np.uint8)
    old12=bits(old_flags[11*B:]);current=raw_groups.get(('broad',12,'all'))
    comparison=None
    if current:comparison=dict(old179_correct=nbits(old12),new_correct=nbits(current['candidate']),
        newly_correct=nbits(current['candidate']&~old12),previously_correct_lost=nbits(old12&~current['candidate']))
    summary=dict(reference_assisted_partition=True,completed_new_configurations=len(new),expected_new_configurations=319,completed_ancillary_raw_configurations=len(ancillary),
        ancillary_scope=json.loads((OUT/'raw_ancillary_scope.json').read_text()) if (OUT/'raw_ancillary_scope.json').exists() else None,
        groups=summary_groups,block12_comparison_to_original_179=comparison,
        design_assistance='Full design sequences are used to assign reads to blocks; reads are not corrected, filled, or converted to reference sequences',
        metrics='Candidate unions and output unions are truth-labelled coverage; each decision output is chosen without truth before scoring',
        new_read_limit=200,live_path_limit=1050,global_index_domain=N,
        fair_partition_control='assigned_pooled versus broad/common10 uses the exact same aligned/filtered read population, joined versus separate',
        split_counts=json.loads((OUT/'reads/split_counts.json').read_text()) if (OUT/'reads/split_counts.json').exists() else None,
        pilot_validation=json.loads((OUT/'pilot_split/sequence_validation.json').read_text()) if (OUT/'pilot_split/sequence_validation.json').exists() else None,
        full_sequence_quality_audit=json.loads((OUT/'full_fastq_audit.json').read_text()) if (OUT/'full_fastq_audit.json').exists() else None,
        assigned_pool_audit=json.loads((OUT/'assigned_pool_audit.json').read_text()) if (OUT/'assigned_pool_audit.json').exists() else None)
    if args.partial:
        (OUT/'partial_summary.json').write_text(json.dumps(summary,indent=2));print('PARTIAL',len(new));return
    # Select one actual decision output for each block. Selection is exploratory
    # using measured lane budgets; the decoder and resulting archive are checked.
    chosen=[raw_groups[('broad',b,'all')]['best_rs'] for b in range(1,13)]
    selected_dir=OUT/'selected_outputs';selected_dir.mkdir(exist_ok=True)
    for row in chosen:
        b=row['block'];selected=selected_dir/f'block{b:02}.strands';fasta=selected_dir/f'block{b:02}.fasta';count=0
        with selected.open('w') as dst,fasta.open('w') as fa,Path(row['path']).open() as src:
            for line in src:
                field,dna=line.rstrip().split('\t');i=int(field)
                if (b-1)*B<=i<b*B:
                    dst.write(line);fa.write(f'>{i}\n{dna}\n');count+=1
        assert count==row['emitted'];row['selected_strands']=str(selected);row['selected_fasta']=str(fasta)
    summary['chosen_per_block']=chosen
    summary['all_block_candidate_union_correct']=sum(nbits(raw_groups[('broad',b,'all')]['candidate']) for b in range(1,13))
    eligible=all(r['within_rs_bound'] for r in chosen)
    summary['file_decode']=dict(attempted=eligible,blocks_in_bound=sum(r['within_rs_bound'] for r in chosen))
    if eligible:
        assembled=OUT/'combined_actual_outputs.strands'
        with assembled.open('w') as handle:
            for b,row in enumerate(chosen,1):
                handle.write(Path(row['selected_strands']).read_text())
        archive=OUT/'file3_reference_assisted_recovered.zip'
        with (OUT/'outer_decode.log').open('w') as log:p=subprocess.run([str(DATA/'bin/rs_driver_fast'),'decode',str(assembled),str(archive),'12'],stdout=log,stderr=log,timeout=600)
        summary['file_decode']['exit_code']=p.returncode
        if not p.returncode:
            content=archive.read_bytes();exact=content==(DATA/'results/reference_file3.zip').read_bytes();crc=zipfile.ZipFile(archive).testzip()
            summary['file_decode'].update(byte_identical=exact,zip_crc_error=crc,sha256=hashlib.sha256(content).hexdigest(),path=str(archive))
            assert exact and crc is None
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    (OUT/'all_decision_metrics.json').write_text(json.dumps(decisions,indent=2))
    write_tsv(OUT/'candidate_comparison.tsv',configs)
    columns=['mode','block','k','cutoff','ratio','policy','emitted','correct','errors','erasures','max_lane_budget','within_rs_bound','precision','recall','outside_target_block_outputs','path']
    with (OUT/'decision_comparison.tsv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(decisions)
    if current:
        with (OUT/'block12_per_index.tsv').open('w') as handle:
            handle.write('index\toriginal179_success\tseparated179_success\tnewly_correct\tlost_previous\n')
            for j in range(B):
                before=bool(old12&(1<<j));after=bool(current['candidate']&(1<<j))
                handle.write(f'{11*B+j}\t{int(before)}\t{int(after)}\t{int(after and not before)}\t{int(before and not after)}\n')
    plot(summary)
    write_report(summary)
    files=['summary.json','REPORT.zh-CN.md','candidate_comparison.tsv','decision_comparison.tsv','block12_per_index.tsv','comparison.png','comparison.svg']
    if (OUT/'full_fastq_audit.json').exists():files.append('full_fastq_audit.json')
    if (OUT/'assigned_pool_audit.json').exists():files.append('assigned_pool_audit.json')
    validation=dict(configurations=319,completed_ancillary_raw_configurations=len(ancillary),pilot_validation=summary['pilot_validation'],full_sequence_quality_audit=summary['full_sequence_quality_audit'],assigned_pool_audit=summary['assigned_pool_audit'],
        read_limit60_control=json.loads((OUT/'read_limit_control.json').read_text()),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'benchmarks/antkowiak2020/split_aligned_blocks.cpp',DATA/'bin/graph_block_rs']},
        artifact_sha256={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in files})
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2));print('BLOCK DECODE REPORT COMPLETE',flush=True)

def write_tsv(path,rows):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)

def plot(s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lookup={(r['mode'],r['block'],r['subset']):r for r in s['groups']}
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained')
    x=np.arange(1,13)
    for j,(mode,label) in enumerate([('historical_processed','Original processed pool'),('assigned_pooled','Assigned reads re-pooled'),('broad','Assigned reads separated')]):
        y=[lookup[(mode,b,'common10')]['candidate_union_correct'] for b in x]
        axes[0].bar(x+(j-1)*.25,y,.25,label=label)
    axes[0].axhline(10977,color='black',ls='--',linewidth=1,label='10,977: erasure-only threshold')
    axes[0].set(xlabel='Block',ylabel='Exact candidate union (same 10 settings)',xticks=x,title='A. Matched-grid reference-assisted experiment')
    axes[0].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.14),ncol=2)
    labels=[];candidate=[];correct=[];errors=[]
    for mode,subset,label in [('strict','common10','Strict / 10'),('inclusive','common10','Inclusive / 10'),('broad','common10','Broad / 10'),('broad','all','Broad / 179')]:
        r=lookup[(mode,12,subset)];labels.append(label);candidate.append(r['candidate_union_correct'])
        correct.append(r['best_single_unique']['correct']);errors.append(r['best_single_unique']['errors'])
    xx=np.arange(4)
    axes[1].bar(xx-.2,candidate,.4,label='Candidate union across settings')
    axes[1].bar(xx+.2,correct,.4,label='Best unique output: correct')
    axes[1].bar(xx+.2,errors,.4,bottom=correct,label='Best unique output: errors',color='#c44e52')
    axes[1].axhline(10977,color='black',ls='--',linewidth=1)
    axes[1].set(xticks=xx,xticklabels=labels,ylabel='Block 12 strands',title='B. Coverage union versus actual unique decisions')
    axes[1].legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.14))
    fig.savefig(OUT/'comparison.png',dpi=180);fig.savefig(OUT/'comparison.svg');plt.close(fig)

def write_report(s):
    lookup={(r['mode'],r['block'],r['subset']):r for r in s['groups']};c=s['split_counts'];d=s['block12_comparison_to_original_179']
    lines=['**File 3：与设计序列比对分块，再独立 DBGPS2＋RS 解码**','',
        '本实验先利用完整设计片段进行 reads 归属，再使用原始 reads 做图搜索及内码纠错。分块使用了已知 payload，因此属于参考辅助实验，不能作为无参考端到端恢复率。比对只决定归组，不修改 reads 碱基、质量值或补齐序列。', '',
        f"全部 {c['total_reads']:,} 条原始 reads 均参与比对；归入 broad 块池 {sum(r['broad_reads'] for r in c['blocks']):,} 条，跨块竞争 {c['ambiguous_reads']:,} 条，未达到归组要求 {c['unassigned_reads']:,} 条。低置信与跨块歧义分别保留。",'',
        '比对：minimap2 k11/w2，允许替换与indel；最多5条次优比对。broad 要求至少20个匹配、30个已比对query碱基、20nt payload重叠、AS≥25、匹配占比≥60%，且得分不低于最佳合格得分80%的已报告合格命中必须属于同一块。strict 在此基础上要求主比对 MAPQ≥10、至少24个匹配、35nt query比对、25nt payload重叠、AS≥40及匹配占比≥70%。未报告的替代比对仍可能存在，因此 broad 不是归属正确性的保证。','',
        '**逐块 reads 与共同10组参数的正确候选并集**','',
        '| 块 | broad reads | strict reads | 旧处理后混池 | 同批已分配reads混池 | 分块独立 |', '|---:|---:|---:|---:|---:|---:|']
    for row in c['blocks']:
        b=row['block'];counts=[lookup[(mode,b,'common10')]['candidate_union_correct'] for mode in ['historical_processed','assigned_pooled','broad']]
        lines.append(f"| {b} | {row['broad_reads']:,} | {row['strict_reads']:,} | "+' | '.join(f'{x:,}' for x in counts)+' |')
    lines+=['','共同参数：k15/c1/r1.8、2.55；k17/18/19各c1/r2.55、3.5；k20/c1/r3.5、5.0。旧处理后混池保留55–70nt并截至60nt；本次原始输入保留实际长度，计数上限200nt。所有实验使用相同1050活跃路径上限、原始全局index与相同内码RS单符号纠错。','',
        f"![分块与重混对照]({OUT/'comparison.png'})",'',
        '图中虚线只表示没有错误输出时的10977条门槛；真实外码可恢复性以下方六列2E+S预算为准。','',
        '“同批已分配reads混池”是把12个broad文件直接拼接重混。它与“分块独立”使用相同reads，能区分分块降噪与比对筛选/长度预处理的作用。该输入一致性另以全记录SHA-256多重集独立核对。候选正确并集按完整60nt真值和纠错后index去重；正确候选存在不等于能无误选出。','',
        f"核心实验完整完成319组新配置。另有未筛选全原始混池的 {s['completed_ancillary_raw_configurations']} 组已完成辅助配置，单列于完整统计和逐参数表；该辅助对照未完成全部10组，部分高噪声搜索主动中止，不将中止配置记为零，不纳入共同10组比较。",'',
        '**第12块额外扫描全部179组参数**','',
        f"原处理后混池179组正确候选并集为 **{d['old179_correct']:,}**；参考辅助broad分块179组为 **{d['new_correct']:,}**。新增正确索引 {d['newly_correct']:,}，此前正确但本次未覆盖 {d['previously_correct_lost']:,}。两次读集与分组前处理不同，属于整体流程比较。",'',
        '| 归组方式 | 参数数 | 正确候选并集 | 最佳单参数正确候选 | 最佳最终单条输出：正确/错误 |', '|---|---:|---:|---:|---|']
    for mode,subset in [('broad','common10'),('strict','common10'),('inclusive','common10'),('broad','all')]:
        r=lookup[(mode,12,subset)];u=r['best_single_unique']
        lines.append(f"| {mode}/{subset} | {r['parameters']} | {r['candidate_union_correct']:,} | {r['best_single_candidate']['correct']:,} | {u['correct']:,} / {u['errors']:,} |")
    old_k_path=DATA/'results/analyzer_block12/summary.json'
    old_k={r['k']:r['per_block'][11]['correct_after_rs_any'] for r in json.loads(old_k_path.read_text())['graph_feasibility_by_k']} if old_k_path.exists() else {}
    lines+=['','**第12块逐k的参数并集**','', '| k | 原处理后混池 | 参考辅助分块、保留原长度 |','|---:|---:|---:|']
    for k in [13,15,17,18,19,20,21,25,31]:
        r=lookup[('broad',12,'k'+str(k))]
        lines.append(f"| {k} | {old_k.get(k,'—')} | {r['candidate_union_correct']:,} |")
    lines+=['','inclusive 还纳入已报告近优比对支持第12块的跨块歧义reads；这些reads仍保留在独立歧义文件中。此变体只扫描共同10组参数。','',
        '**每块实际唯一输出与外码预算**','',
        '| 块 | k/c/r/规则 | 正确 | 错误 | 缺失 | 最大列2E+S | 满足5406界 |','|---:|---|---:|---:|---:|---:|---|']
    for r in s['chosen_per_block']:
        lines.append(f"| {r['block']} | {r['k']}/{r['cutoff']}/{r['ratio']:.2f}/{r['policy']} | {r['correct']:,} | {r['errors']:,} | {r['erasures']:,} | {r['max_lane_budget']:,} | {'是' if r['within_rs_bound'] else '否'} |")
    lines+=['','每块选择已测试输出中最大列预算最小的一组参数/规则，选择利用真值评估，属于探索性参数选择；实际输出的候选选择本身不读真值。外码预算逐6个GF(2^14)码字计算。', '']
    f=s['file_decode']
    if f.get('byte_identical'):lines += [f"完整文件通过实际外码解码、ZIP CRC及与参考ZIP逐字节一致校验：[恢复文件]({f['path']})。这是参考辅助分块条件下的文件恢复。"]
    else:lines += [f"{f['blocks_in_bound']}/12块存在满足外码纠错界的已测实际输出；未获得通过验证的新增完整File3。"]
    lines+=['','**文件**','',
        f"- [第12块broad FASTQ]({OUT/'reads/broad/block12.fq.gz'})；[strict FASTQ]({OUT/'reads/strict/block12.fq.gz'})；[含歧义支持的FASTQ]({OUT/'reads/block12.inclusive.fq.gz'})。",
        f"- [第12块选定配置的实际片段输出]({OUT/'selected_outputs/block12.strands'})；[FASTA]({OUT/'selected_outputs/block12.fasta'})。该输出只按块地址过滤，保留其中的错误链；没有用真值删除错误。其余块同名编号01至12。",
        f"- 其余块文件在 `{OUT/'reads/broad'}` 和 `{OUT/'reads/strict'}`，命名 `block01.fq.gz` 至 `block12.fq.gz`。",
        f"- [第12块逐索引对照]({OUT/'block12_per_index.tsv'})；[逐参数候选结果]({OUT/'candidate_comparison.tsv'})；[所有唯一输出指标]({OUT/'decision_comparison.tsv'})。",
        f"- [分块规则、完整命令和计数]({OUT/'split_provenance.json'})；[完整统计]({OUT/'summary.json'})；[校验]({OUT/'validation.json'})。",'',
        f"验证：pilot检查记录为 `{s['pilot_validation']}`；全量reads计数守恒；新实验程序在read limit=60时的11种决策输出均与原程序逐字节一致。生产DBGPS2及原实验二进制未改动。",
        '',f"全量原始序列和质量多重集核对：{'通过' if s['full_sequence_quality_audit'] and s['full_sequence_quality_audit'].get('passed') else '结果尚未生成'}。[独立核对记录]({OUT/'full_fastq_audit.json'})。"]
    (OUT/'REPORT.zh-CN.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
