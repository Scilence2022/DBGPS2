#!/usr/bin/env python3
"""Collect completed experiments, verify recoveries, and write a Chinese report."""
import csv
import hashlib
import json
import pathlib
import subprocess
import sys
import time
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT/'input/antkowiak2020'
OUT = DATA/'results'

def load(path):
    return json.loads(path.read_text())

def truth(label):
    p=DATA/f'sources/File{label[-1]}_design.fasta' if label=='file1' else DATA/f'sources/File{label[-1]}_synthesized.fasta'
    return [s for s in p.read_text().splitlines() if s and not s.startswith('>')]

def score(path, target):
    outputs={int(x.split('\t')[0]):x.split('\t')[1] for x in path.read_text().splitlines()}
    correct=sum(i<len(target) and s==target[i] for i,s in outputs.items())
    blocks=[]
    for start in range(0,len(target),16383):
        end=min(len(target),start+16383)
        emitted=sum(start<=i<end for i in outputs)
        errors=sum(start<=i<end and s!=target[i] for i,s in outputs.items())
        blocks.append(dict(block=start//16383,errors=errors,erasures=end-start-emitted,
                           budget=2*errors+end-start-emitted,within_rs_bound=2*errors+end-start-emitted<=5406))
    return dict(emitted=len(outputs),correct=correct,errors=len(outputs)-correct,
                erasures=len(target)-len(outputs),recall=correct/len(target),per_block=blocks)

def verify_decode(label, candidates, expected, blocks, prefix, metrics):
    result=OUT/(prefix+'_decode.json')
    if not all(b['within_rs_bound'] for b in metrics['per_block']):
        record=dict(attempted=False,reason='one or more blocks exceed the conventional RS errors/erasures bound')
    else:
        target=OUT/(prefix+'.zip')
        start=time.monotonic()
        with (OUT/(prefix+'_decode.log')).open('w') as log:
            p=subprocess.run([str(DATA/'bin/rs_driver_fast'),'decode',str(candidates),str(target),str(blocks)],stdout=log,stderr=log,timeout=600)
        record=dict(attempted=True,exit_code=p.returncode,seconds=time.monotonic()-start)
        if p.returncode==0:
            b=target.read_bytes()
            record.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),
                          byte_identical=b==expected.read_bytes(),zip_error=zipfile.ZipFile(target).testzip())
            assert record['byte_identical'] and record['zip_error'] is None, record
    result.write_text(json.dumps(record,indent=2))
    return record

def main():
    required=[OUT/f'{label}_{kind}.json' for label in ['file1','file2','file3'] for kind in ['graph_metrics','msa']]
    required.append(OUT/'file1_notebook.json')
    while not all(p.exists() and p.stat().st_size for p in required):
        time.sleep(10)
    all_rows={label:load(OUT/f'{label}_graph_metrics.json') for label in ['file1','file2','file3']}
    selected={label:max(rows,key=lambda r:(r['correct'],-r['errors'],-r['overflow'])) for label,rows in all_rows.items()}
    summary={}
    for label,best in selected.items():
        target=truth(label)
        eligible=[r for r in all_rows[label] if all(b['within_rs_bound'] for b in r['per_block'])]
        file_choice=max(eligible,key=lambda r:(r['correct'],-r['errors'])) if eligible else best
        name=f"{label}.k{file_choice['k']}.c{file_choice['cutoff']}.r{file_choice['ratio']:.2f}.strands"
        expected=DATA/'sources/noisy_dna_data_storage/data/File_1.zip' if label=='file1' else OUT/f'reference_{label}.zip'
        decoded=verify_decode(label,OUT/name,expected,len(target)//16383,label+'_best',file_choice)
        decoded['parameters']={k:file_choice[k] for k in ['k','cutoff','ratio']}
        prefix=OUT/(label+'_msa_inner')
        with prefix.with_suffix('.json').open('w') as log:
            subprocess.run([str(DATA/'bin/inner_select'),str(OUT/(label+'_msa.candidates.txt')),str(len(target)),str(prefix)],stdout=log,check=True)
        inner={mode:score(OUT/(label+'_msa_inner.'+mode+'.strands'),target) for mode in ['strict','corrected']}
        summary[label]=dict(graph=best,msa=load(OUT/(label+'_msa.json')),msa_inner=inner,
                            reads=load(OUT/(label+'_reads.json')),decode=decoded,
                            best_block_feasibility=max(all_rows[label],key=lambda r:(sum(b['within_rs_bound'] for b in r['per_block']),r['correct'],-r['errors'])),
                            configurations_with_all_blocks_in_bound=len(eligible))
    # The released notebook is a separate preprocessing/parameter baseline.
    notebook=load(OUT/'file1_notebook.json')
    with (OUT/'file1_notebook_inner.json').open('w') as log:
        subprocess.run([str(DATA/'bin/inner_select'),str(OUT/'file1_notebook.candidates.txt'),'16383',str(OUT/'file1_notebook_inner')],stdout=log,check=True)
    notebook_inner={mode:score(OUT/('file1_notebook_inner.'+mode+'.strands'),truth('file1')) for mode in ['strict','corrected']}
    summary['file1']['notebook']=notebook
    summary['file1']['notebook_inner']=notebook_inner
    # A matched postprocessing comparison: both candidate generators use the
    # same GF(64) inner decoder and error/erasure budget accounting.
    nmetrics=notebook_inner['strict']
    summary['file1']['notebook_decode']=verify_decode('file1',OUT/'file1_notebook_inner.strict.strands',DATA/'sources/noisy_dna_data_storage/data/File_1.zip',1,'file1_notebook',nmetrics)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    flat=[r for rows in all_rows.values() for r in rows]
    columns=['dataset','k','cutoff','ratio','truth','emitted','correct','errors','erasures','exact_recall','precision','overflow','ambiguous','graph_supported']
    with (OUT/'graph_comparison.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(flat)
    plot(all_rows,summary)
    write_report(all_rows,summary)
    print('REPORT READY',OUT/'REPORT.zh-CN.md',flush=True)

def plot(rows,summary):
    sys.path.insert(0,str(DATA/'plotlibs'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    colors=['#176b97','#ed9b40','#42966d']
    for label,color in zip(['file1','file2','file3'],colors):
        ks=sorted(set(r['k'] for r in rows[label]));best=[max((r for r in rows[label] if r['k']==k),key=lambda r:r['correct']) for k in ks]
        axes[0].plot(ks,[r['exact_recall']*100 for r in best],'-o',label=label.title(),color=color,linewidth=2)
    axes[0].set(xlabel='k-mer size',ylabel='Exact strand recall (%)',title='Exploratory parameter sweep',ylim=(-2,103))
    axes[0].legend(frameon=False);axes[0].grid(alpha=.2)
    axes[0].set_xticks(sorted({r['k'] for data in rows.values() for r in data}))
    axes[0].tick_params(axis='x',labelsize=8)
    labels=['File 1','File 2','File 3'];x=np.arange(3);w=.24
    for offset,key,color,title in [(-1,'raw','#a3adb8','Observed'),(0,'msa','#ed9b40','MSA + vote'),(1,'graph','#176b97','Adapted DBGPS2')]:
        values=[]
        for label in ['file1','file2','file3']:
            d=summary[label]
            values.append(100*(d['msa']['raw_exact_recovered']/d['graph']['truth'] if key=='raw' else d['msa']['recall'] if key=='msa' else d['graph']['exact_recall']))
        bars=axes[1].bar(x+offset*w,values,w,color=color,label=title)
        axes[1].bar_label(bars,fmt='%.1f',fontsize=8,padding=3)
    axes[1].set(xticks=x,xticklabels=labels,ylim=(0,109),ylabel='Exact strand coverage (%)',title='Same retained sequencing reads')
    axes[1].legend(frameon=False,fontsize=7,loc='upper center',bbox_to_anchor=(.5,-.11),ncol=3)
    file1=rows['file1'];ks=sorted(set(r['k'] for r in file1));best=[max((r for r in file1 if r['k']==k),key=lambda r:r['correct']) for k in ks]
    axes[2].plot(ks,[100*r['graph_supported']/r['truth'] for r in best],'-s',color='#8b6bb1',label='All true k-mers exceed cutoff')
    axes[2].plot(ks,[100*r['exact_recall'] for r in best],'-o',color='#176b97',label='Unique, correct RS-valid output')
    axes[2].set(xlabel='k-mer size',ylabel='Fraction of File 1 strands (%)',ylim=(-2,103),title='Connectivity and ambiguity')
    axes[2].legend(frameon=False,fontsize=8);axes[2].grid(alpha=.2)
    axes[2].set_xticks([12,13,15,17,19,21,25,31])
    fig.savefig(OUT/'comparison.png',dpi=180)
    fig.savefig(OUT/'comparison.svg')
    plt.close(fig)

def write_report(rows,summary):
    pct=lambda x:f'{100*x:.2f}%'
    text=['**Antkowiak et al. (2020)：编码分析与 DBGPS2 片段恢复实测**','',
          '本报告使用完整公开测序数据，所有原始 FASTQ 均通过公布的 MD5。结果中的“DBGPS2”均指实验适配程序：复用本仓库的精确 k-mer 计数、规范化图查询及动态覆盖比规则，增加无引物的中部索引双向搜索和论文所用 GF(64) 校验。原版 CLI 不能直接解释本数据格式。','',
          '**编码与适配**','',
          '每条链包含 42 nt 数据（84 bit）、12 nt 索引（24 bit）和 6 nt 内层 RS 校验（12 bit），共 60 nt。外层 RS 为 (16383,10977)/GF(2^14)，每条链承载 6 个外码符号；内层为 (20,18)/GF(2^6)。碱基映射为 A=00、C=01、G=10、T=11。公开实现采用有限域加法进行伪随机化；这批实测数据的索引掩码为 1071068。完整文件仍需要外层纠错。[论文](https://www.nature.com/articles/s41467-020-19148-3) [作者代码](https://github.com/MLI-lab/noisy_dna_data_storage)','',
          'DBGPS2 的默认前端索引、固定引物和 CRC16/字节 RS 与该格式不同；42 nt 数据区的 84 bit 也不满足原版参数要求的整字节布局。因此不能只换一个 `-M` 或 `-E rs` 参数。适配程序仅从已知索引和读段图生成候选，不把参考数据区当作种子。每代路径上限保持 1050；溢出作为缺失处理。严格 GF(64) 校验后仅接受唯一结果。','',
          '**数据与参考核验**','',
          '| 数据集 | 原始 reads | 保留 reads（55–70 nt，截至 60 nt） | 设计链数 |','|---|---:|---:|---:|']
    for label in ['file1','file2','file3']:
        d=summary[label];text.append(f"| {label.title()} | {d['reads']['total_reads']:,} | {d['reads']['filtered_reads']:,} | {d['graph']['truth']:,} |")
    text += ['', 'File 1 的原始设计序列及 ZIP 来自作者仓库。File 2/3 的设计文件来自后续公开分析仓库：它们把索引移到前端并重新生成校验，不能直接与物理实验的完整读段比较。本测试保留其公开的 42 nt 数据区，按实测编码重排索引、加入原随机化并重算内层校验，生成评估参考。该转换精确重现 File 1 全部 16383 条真值；File 2/3 转换后的参考可无错误通过外层解码并生成内部 CRC 完整的 ZIP，且在原始 reads 中找到完整精确匹配。参考转换只用于评估，参考文件由恢复程序在搜索结束后读取。[后续数据来源](https://github.com/fml-ethz/dt4dds-challenges_notebooks)','',
             '**完整数据的片段比较**','',
             '以下使用与论文相同的“至少得到一次完全正确的 60 nt 片段”定义。图方法还列出错误输出数。MSA 基线使用相同保留 reads，prefix=16、簇大小至少 5、MUSCLE 3.8.31、每次最多 15 条、较大簇取 5 组、gap 权重 0.4；随机种子为 20261017。','',
             '| 数据集 | 直接观察到的精确片段 | MSA＋加权投票 | 适配 DBGPS2 | k / 截断 / ratio | 错误输出 |','|---|---:|---:|---:|---|---:|']
    for label in ['file1','file2','file3']:
        d=summary[label];g=d['graph'];text.append(f"| {label.title()} | {pct(d['msa']['raw_exact_recovered']/g['truth'])} | {pct(d['msa']['recall'])} | {g['correct']:,}/{g['truth']:,}（{pct(g['exact_recall'])}） | {g['k']} / {g['cutoff']} / {g['ratio']:.2f} | {g['errors']:,} |")
    n=summary['file1']['notebook'];text += ['',f"另按作者发布 notebook 的设置（prefix=14、最小簇 4、gap 权重 0.5，先对 55–70 nt 完整 reads 比对、投票后截取 60 nt）重跑 File 1，准确覆盖 {n['recovered']:,}/16383，即 {pct(n['recall'])}。它与同输入基线的预处理不同。论文报告的 naive/LSH 数值分别为约 74%/83%（File 1），LSH 对 File 2/3 约 35%/31%；这些是文献值，本次没有重跑完整 LSH。",'',
             '为区分候选覆盖和最终片段决策，又对 MSA 候选统一执行 GF(64) 内码和按索引投票。`strict` 只保留无需修正的码字；`corrected` 允许修正一个 6-bit 符号，并优先选择修正次数更少的候选；票数并列时视为缺失。','',
             '| 数据集 | MSA strict：正确 / 错误 | MSA corrected：正确 / 错误 |','|---|---:|---:|']
    for label in ['file1','file2','file3']:
        a=summary[label]['msa_inner'];text.append(f"| {label.title()} | {a['strict']['correct']:,} / {a['strict']['errors']:,} | {a['corrected']['correct']:,} / {a['corrected']['errors']:,} |")
    text += ['',f"![恢复比较]({OUT/'comparison.png'})",'',
             '**k 与覆盖阈值**','',
             '较小 k 提高正确局部片段的存活率，但也把许多错误及不同原链混合到同一图中。因此小 k 必须与覆盖截断共同调节。极小 k 下，即使真值路径在图中存在，也可能因分支超过上限而无法唯一恢复；较大 k 则丢失更多低覆盖的正确路径。图方法可把不同错误读段中仍正确的局部片段组合起来，并在路径搜索中利用索引和内码筛选；MSA 基线依赖分簇、对齐和投票产生正确的全长候选。','',
             '| File 1 的 k | 本扫描中最佳截断 | 精确恢复率 | 错误输出 | 路径溢出 |','|---:|---:|---:|---:|---:|']
    for k in sorted(set(r['k'] for r in rows['file1'])):
        r=max((r for r in rows['file1'] if r['k']==k),key=lambda r:r['correct']);text.append(f"| {k} | {r['cutoff']} | {pct(r['exact_recall'])} | {r['errors']} | {r['overflow']} |")
    text += ['', '表中最佳参数根据本数据的参考真值选出，属于探索性参数扫描；没有独立留出集，不能据此保证同参数在其他实验中最优。每次扫描的完整参数及逐索引状态均保留。','',
             '**最终文件验证**','',
             '每个外码块的常规错误/缺失纠错条件为 `2E + S <= 5406`。只有所有块满足条件后才尝试完整解码；最终还要求字节一致及 ZIP 内部 CRC 通过。外层使用作者的频域消息解码方式和旧版随机化；GF(2^14) 乘法、求逆和幂采用等价的对数表加速，已与原实现进行 232072 项代数校验，并用 File 1 原 ZIP 做字节级交叉验证。','',
             '| 数据集 | 满足外码界的块 | 测序恢复文件 |','|---|---:|---|']
    for label in ['file1','file2','file3']:
        d=summary[label];b=d['graph']['per_block'];decoded=d['decode'];result='通过字节及 ZIP 校验' if decoded.get('byte_identical') else '未取得通过验证的完整文件'
        text.append(f"| {label.title()} | {sum(x['within_rs_bound'] for x in b)}/{len(b)} | {result} |")
    last=summary['file3'];worst=max(last['graph']['per_block'],key=lambda b:b['budget'])
    if not last['decode'].get('byte_identical'):
        text += ['',f"File 3 的整体精确恢复率为 {pct(last['graph']['exact_recall'])}，但外码逐块工作。最弱的第 {worst['block']+1} 块缺失 {worst['erasures']:,}/16383 条（{pct(worst['erasures']/16383)}），另有 {worst['errors']} 条错误，`2E+S={worst['budget']}` 超过 5406。当前输出对该块没有足够片段完成原外码恢复；整体平均比例不能代替逐块判断。"]
    block_best=last.get('best_block_feasibility')
    if block_best and sum(b['within_rs_bound'] for b in block_best['per_block'])>sum(b['within_rs_bound'] for b in last['graph']['per_block']):
        text += ['',f"另按达标块数量选择参数，File 3 最多有 {sum(b['within_rs_bound'] for b in block_best['per_block'])}/12 块满足外码界，对应 k={block_best['k']}、c={block_best['cutoff']}、ratio={block_best['ratio']:.2f}，精确恢复率 {pct(block_best['exact_recall'])}。它与整体恢复率最高的配置不同，仍不足以恢复完整文件。"]
    text += ['', 'File 1 与作者公开原 ZIP 比较；File 2/3 与从完整公开参考序列解出的文件比较，且后者也通过 ZIP 内部 CRC。参考文件解出的 ZIP 与测序恢复的 ZIP 具有不同文件名，避免混淆。','']
    for label in ['file1','file2','file3']:
        d=summary[label]['decode']
        if d.get('byte_identical'):text.append(f"{label.title()}：{d['bytes']:,} 字节，SHA-256 `{d['sha256']}`。")
    text += ['', '**复现与边界**','',
             f"完整指标：[summary.json]({OUT/'summary.json'})；全部图参数：[graph_comparison.csv]({OUT/'graph_comparison.csv'})；参考格式转换：[reference_transcoding.json]({OUT/'reference_transcoding.json'})；可复现程序：[README]({ROOT/'benchmarks/antkowiak2020/README.md'})。原始数据保存在 `{DATA/'raw'}`，参考及上游代码保存在 `{DATA/'sources'}`。",'',
             '本结果对应本仓库版本、完整 reads、上述预处理及有限参数网格。恢复与评估分离；原链数据未参与图搜索。MSA 仅运行一个固定随机种子，未计算重复实验置信区间。运行时间受同时进行的扫描与多进程负载影响，仅作为本机记录，不作为跨硬件速度结论。未修改生产版 DBGPS2 源文件，未提交或推送任何内容。']
    blind=OUT/'reference_blind_check.json'
    if blind.exists() and load(blind).get('output_identical'):
        text += ['',f"独立检查：对 File 1 完全不传入参考 FASTA，仍得到相同的 16125 条输出，字节级结果一致。[检查记录]({blind})。"]
    supplement=OUT/'FILE3_RATIO_REPORT.zh-CN.md'
    if supplement.exists():
        text += ['',f"File 3 的小 k / ratio 联合扫描及图支持上限分析：[补充报告]({supplement})。"]
    (OUT/'REPORT.zh-CN.md').write_text('\n'.join(text)+'\n')

if __name__=='__main__':
    main()
