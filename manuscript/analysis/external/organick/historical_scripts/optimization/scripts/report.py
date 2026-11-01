#!/usr/bin/env python3
"""Create the full-read optimization report from completed, measured outputs."""
import csv
import gzip
import json
from pathlib import Path
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2];OPT=ROOT/'optimization'
ensemble=json.loads((OPT/'results/refined.header13.ensemble_summary.json').read_text())
bound=json.loads((OPT/'results/graph_bound_union.json').read_text())['bounds'][0]['either_k_either_direction']
native=json.loads((OPT/'results/full_summary.json').read_text())
inner=json.loads((OPT/'results/inner_full_summary.json').read_text())
truth=set()
with gzip.open(ROOT/'data/id20.refs.txt.gz','rt') as handle:
    for line in handle:truth.add(line.strip()[3:-3][20:130])
nohp=re.compile('AA|CC|GG|TT')
profiles=[]
for kind,rows in [('native_search',native),('inner_validated',inner)]:
    for row in rows:
        predictions={}
        with open(row['output']) as handle:
            for candidate in csv.DictReader(handle,delimiter='\t'):
                body=candidate['payload']
                if nohp.search(body):continue
                key=body[:13];old=predictions.get(key)
                if old is None:predictions[key]=body
                elif old!=body:predictions[key]='CONFLICT'
        accepted=[body for body in predictions.values() if body!='CONFLICT']
        exact=sum(body in truth for body in accepted)
        profiles.append({'kind':kind,'tag':row['tag'],'exact':exact,'returned':len(accepted),
            'wrong':len(accepted)-exact,'exact_rate':exact/607150,'precision':exact/len(accepted),
            'k':row['k'],'cutoff':row['cutoff'],'ratio':row['ratio'],'skip':row['skip']})
best_native=max((x for x in profiles if x['kind']=='native_search'),key=lambda row:(row['exact'],-row['wrong']))
best_inner=max((x for x in profiles if x['kind']=='inner_validated'),key=lambda row:(row['exact'],-row['wrong']))
maximum=ensemble['selected']['maximum_exact'];balanced=ensemble['selected']['precision_99_99'];zero=ensemble['selected']['zero_observed_errors']
assert maximum['exact']<=bound
for name,row in ensemble['selected'].items():
    if row is None:continue
    sequences=[]
    with gzip.open(row['fasta'],'rt') as handle:
        for line in handle:
            if not line.startswith('>'):sequences.append(line.strip())
    assert len(sequences)==row['returned'] and len(set(s[:13] for s in sequences))==len(sequences)
    assert all(len(s)==110 for s in sequences)
    assert sum(s in truth for s in sequences)==row['exact']
    assert sum(s not in truth for s in sequences)==row['wrong']
print('All selected FASTA sequence counts, full-length exactness and uniqueness verified.',flush=True)
baseline={'exact':594630,'wrong':25,'returned':594655,'exact_rate':594630/607150,'precision':594630/594655}
labels=['前次全量基线','仅调参＋双向','内码校验单配置','组合：最多精确恢复','组合：正确率≥99.99%','组合：观察到零错误']
display=[baseline,best_native,best_inner,maximum,balanced,zero]
table=['| 策略 | 输出条数 | 精确恢复条数 | 精确恢复率 | 错误输出 | 输出正确率 |','|---|---:|---:|---:|---:|---:|']
for label,row in zip(labels,display):
    table.append(f'| {label} | {row["returned"]:,} | {row["exact"]:,} | {row["exact_rate"]:.5%} | {row["wrong"]:,} | {row["precision"]:.5%} |')
metrics=[dict(strategy=label,**row) for label,row in zip(labels,display)]
(OPT/'results/final_summary.json').write_text(json.dumps({'rows':metrics,'path_existence_bound':bound,'rs_validation':'not_performed'},indent=2))
with (OPT/'results/final_summary.tsv').open('w') as handle:
    fields=['strategy','exact','returned','wrong','exact_rate','precision']
    writer=csv.DictWriter(handle,delimiter='\t',fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(metrics)

fig,axes=plt.subplots(1,2,figsize=(10.5,4.8),layout='constrained')
short=['Baseline','Tuned','Inner code','Max exact','>=99.99%','Zero observed']
colors=['#8B9298','#548AB3','#298F87','#D98435','#5874B0','#7B6AA5']
axes[0].scatter(range(len(short)),[row['exact_rate']*100 for row in display],color=colors,s=65,zorder=3)
axes[0].set_xticks(range(len(short)),short)
for i,row in enumerate(display):
    axes[0].annotate(f'{row["exact_rate"]*100:.3f}%',(i,row['exact_rate']*100),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
axes[0].axhline(bound/607150*100,ls='--',color='#A54848',label=f'Path-existence bound: {bound/607150:.4%}')
axes[0].set_ylim(90,99)
axes[0].set_ylabel('Exact 110 nt recovery (% of 607,150)')
axes[0].set_title('100% of raw reads: 25.8913x coverage')
axes[0].legend(fontsize=8,loc='lower left')
axes[1].bar(short,[row['wrong'] for row in display],color=colors)
axes[1].set_ylabel('Incorrect output count')
axes[1].set_title('Recovery and correctness are separate objectives')
for ax in axes:
    ax.tick_params(axis='x',rotation=30)
    ax.spines[['right','top']].set_visible(False)
    ax.set_axisbelow(True);ax.grid(axis='y',color='#E6E6E6')
fig.savefig(OPT/'results/optimization_comparison.png',dpi=180)
fig.savefig(OPT/'results/optimization_comparison.svg')

def rule_text(row):
    return f'最少 {row["minimum_votes"]} 票，候选一致性至少 {row["minimum_agreement"]:.0%}，领先票差至少 {row["minimum_margin"]}，种子观测支持度至少 {row["minimum_seed_support"]}'
base_meta=json.loads((OPT/'results/provenance.json').read_text())
inner_meta=json.loads((OPT/'results/inner_provenance.json').read_text())
refined=json.loads((OPT/'results/refined_summary.json').read_text())
text=f'''# 原始 reads 使用率 100%：DBGPS2 恢复优化

本次“100%”指使用全部 **15,719,903 条原始 reads**，对应平均 **25.8913×** 原始覆盖；不是 100×，也不是每条目标片段都被无误测到了。设计片段共 **607,150 条**，每条评分区间为完整 110 nt 编码区。

## 结果

在本次已完成的参数和组合规则搜索中，精确恢复条数最大为 **{maximum['exact']:,} 条（{maximum['exact_rate']:.5%}）**，同时输出 **{maximum['wrong']:,} 条错误片段**，输出正确率为 **{maximum['precision']:.5%}**。

如果要求本数据集上观察到零错误，可输出 **{zero['exact']:,} 条（恢复率 {zero['exact_rate']:.5%}）**。这表示已输出片段与设计序列逐碱基相符，不表示全部 607,150 条目标都被恢复，也不是对其他数据的零错误保证。

{chr(10).join(table)}

![全量恢复优化比较]({OPT}/results/optimization_comparison.png)

前次基线为 k=31、cutoff=1、ratio=2.55、skip=0、单向搜索及无同聚物过滤。当前单配置比较使用双向输出，并按 13 nt 编码头拒绝冲突；组合结果按相同编码头投票，不能用参考序列决定某一目标该选择哪条候选。

## 搜索范围与最佳设置

- 全部 runs 均用全部 reads 建图，没有对 reads 降采样。
- 原始候选校验模式：筛选 **{base_meta['screened_configs']}** 个配置。k=21/23/25/27/29/31，cutoff 初始为 0/1/2/3/4/6，ratio 初始为 1.1/1.25/1.4/1.6/1.8/2.0/2.55/3.0，再对领先点局部细化 ratio 和 skip=4/16。
- 为控制搜索成本，初筛查询的是独立固定哈希选择的 **1/128 观测种子**；这些图仍来自全部 reads。挑出的 **{base_meta['full_confirmed_configs']}** 个配置随后对全部 **3,063,745 个正反向种子**实际确认，未将抽样数字外推为最终恢复率。
- 论文内码候选校验模式：另筛选 **{inner_meta['screened_configs']}** 个配置，并全量确认 **{inner_meta['full_confirmed_configs']}** 个配置。
- 再对预测不稳定或尚无稳定预测的正向种子完成 **{len(refined)}** 个细化配置。查询选择只读取已有预测及投票，不按参考真值选择“缺失目标”。
- 原始搜索的最佳全量配置：`{best_native['tag']}`，精确恢复 {best_native['exact']:,} 条、错误 {best_native['wrong']:,} 条。
- 内码校验的最佳单配置：`{best_inner['tag']}`，精确恢复 {best_inner['exact']:,} 条、错误 {best_inner['wrong']:,} 条。

这里的 cutoff 是 DBGPS2 `-c`：cutoff=0 要求 k-mer 计数≥1，cutoff=1 要求计数≥2。

组合使用各配置实际生成的候选，按编码头选择票数最高且不平票的候选，再应用全局一致的票数、领先幅度和种子支持度阈值：

- 最大精确恢复：{rule_text(maximum)}。
- 输出正确率至少 99.99%：{rule_text(balanced)}。
- 观察到零错误：{rule_text(zero)}。

本次共评价 **{ensemble['rules_evaluated']}** 条全局组合规则。不同参数及方向的投票相互关联，票数不应解释成独立证据的概率。

## 格式适配与参考序列的使用

生产 `DBGPS2.c`、建图和覆盖约束路径搜索均未修改。独立实验驱动复用其函数，按实际引物产生种子，并关闭不适用的片段内 CRC/RS。

新增的实验内码验证器放在原有候选验证调用位置：按 `AGCT` 旋转顺序检查相邻碱基转移，27 个四三进制位分组均须小于 64，末尾补齐转移须为 0。它是这套数据的格式适配，不是原生 CLI 的现成开关，也不是 RS 纠错。

旋转顺序和分组相位通过 512 条设计序列做全局标定，随后全部 **607,150 条设计序列均通过该谓词**。前 13 nt 编码头在全部设计序列中也互不重复；组合按这一编码头去重。这里未声称已完整解析论文的二进制 block/column 地址或 RS 布局。

**参考文件用于全局格式标定、参数/组合规则选择，以及最终逐碱基评分；不用于给 reads 分组、提供逐片段起始种子、修补候选或按逐片段答案选路径。** 因为在同一数据集上选优，这属于数据集特定的优化结果，不是独立留出测试，也不是对全部可能参数的数学最优证明。

## 可达路径上界

另做了明确使用真值的诊断：只问“设计序列的正确路径是否完整存在于 reads 建成的图中”，不生成恢复结果。对 k=29/31、现有种子规则、计数至少为 1、取消覆盖比与路径唯一性限制，两个图和方向的并集仅有 **{bound:,} 条（{bound/607150:.5%}）**具备完整路径。

实际最大精确恢复结果距离这个受限的路径存在上界还差 **{bound-maximum['exact']:,} 条**，达到该上界的 **{maximum['exact']/bound:.5%}**。至少 **{607150-bound:,} 条设计片段**在这两个图及现有种子规则下缺少完整可达路径；提高覆盖比或降低 cutoff 不能凭空补出缺失 k-mer。

该数字是有真值指路的图存在性上界，**不能当作 DBGPS2 已经恢复的数量**，也不是其他 k、种子方案、纠错算法或 RS 后的全局上限。可达路径检查已通过无误输入和系统性缺口输入对照。

## RS 验证状态

**本轮仍未执行论文跨片段 RS 纠错。** 所有上述数字都是 RS 前结果。内码谓词和多参数组合不能代替 RS；尚不能声称数字文件已无误恢复。下一层验证需要实际解析内码字节及地址、确定论文的分块和有限域/RS 参数，逐码字纠错后再与设计内容或原文件比较。

## 产物与复现

- `results/refined.header13.maximum_exact.fa.gz`：最多精确恢复策略的实际预测，包含表中残余错误。
- `results/refined.header13.precision_99_99.fa.gz`：输出正确率≥99.99%策略。
- `results/refined.header13.zero_observed_errors.fa.gz`：本数据集上观察到零错误策略。
- `results/final_summary.tsv` / `.json`：最终比较。
- `results/*summary.json`、`*.jsonl`、`*.log`：全部候选配置、筛选、全量确认及细化记录。
- `results/graph_bound_union.json`：独立的图路径上界诊断。
- `results/inner_code_validation.json`、`inner_controls.json`、`bound_controls.json`：格式标定及对照验证。

在仓库根目录执行：

```sh
python3 experiments/nbt4079/optimization/scripts/prepare_seeds.py
uv run --with numpy python experiments/nbt4079/optimization/scripts/validate_inner.py
python3 experiments/nbt4079/optimization/scripts/optimize.py
python3 experiments/nbt4079/optimization/scripts/optimize_inner.py
python3 experiments/nbt4079/optimization/scripts/refine.py
uv run --with numpy python experiments/nbt4079/optimization/scripts/aggregate.py --key-length 13 --include-refinement
python3 experiments/nbt4079/optimization/scripts/compute_bounds.py
uv run --with numpy --with matplotlib python experiments/nbt4079/optimization/scripts/report.py
```

预测 FASTA 已重新读取检查 gzip 完整性、110 nt 长度、编码头唯一性及正确/错误数量。没有提交、推送或合并任何更改。
'''
(OPT/'REPORT.zh-CN.md').write_text(text)
print(OPT/'REPORT.zh-CN.md')
