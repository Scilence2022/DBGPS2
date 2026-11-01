"""Write a candid quantitative report from the completed, locked validation."""
from pathlib import Path
import json
import pandas as pd

R=Path(__file__).resolve().parent


def main():
    p=pd.concat([pd.read_csv(R/'policy_results.csv'),pd.read_csv(R/'crossing_policy_results.csv')],ignore_index=True)
    g=pd.concat([pd.read_csv(R/'native_grid_results.csv'),pd.read_csv(R/'crossing_native_grid_results.csv')],ignore_index=True)
    v=json.loads((R/'validation_summary.json').read_text());tr=pd.read_csv(R/'false_candidate_donor_traces.csv')
    sim=p[p.cohort!='real_HT176_raw'];s=sim[sim.cohort=='scale'];real=p[p.cohort=='real_HT176_raw']
    scale=s.groupby(['k','n','policy'])[['correct_pct','false_accepted','coverage_ceiling_pct']].agg(['mean','std'])
    scale.to_csv(R/'scale_summary.csv')
    error=p[(p.cohort=='error_structure')&(p.policy=='default')].groupby(['k','depth','channel']).correct_pct.agg(['mean','std'])
    error.to_csv(R/'error_structure_summary.csv')
    names={'default':'固定 r=2.55','global_crc_r':'全局 CRC 数量选 r','global_crc_cr':'全局 CRC 数量选 c/r',
           'union_consensus':'一致候选并集','two_support_consensus':'至少两个 r 支持'}
    cv=json.loads((R/'crossing_validation_summary.json').read_text())
    lines=['# DBGPS2 重建边界与无参考选参的深入验证','',
        '**本轮累计完成 366 个模拟计数图、8,784 个参数点和 24 个真实 HT176 参数点。** '
        '其中330个图属于先锁定的主验证，另36个图来自首次跌破阈值边界推导后另行锁定的新文库确认。'
        '新增五种固定决策规则，使用全新随机文库；规模测试最高65,536条设计。第一轮的1,368个图保持不变。'
        '这些是机制与算法行为实验，没有新增速度比较或完整文件解码实验。','',
        '## 结论与证据边界','',
        '当前最稳妥的结论是：**有效连续k窗口的支持概率决定可用的路径证据；局部比例窗口决定这些证据能否在受限搜索中被利用；CRC和参数间一致性仍不能保证所有输出正确。** '
        '模拟把这些环节分开了，但尚未建立经独立实验文库验证的普适损伤阈值。','',
        '## 1 无参考选参的独立验证','',
        '五个规则在新数据测量前冻结。选参模块只能读取每个参数点是否唯一通过原生有效性检查，以及完整候选序列，不能读取设计序列、正确性标签、真路径计数、覆盖上限或模拟参数。'
        '更换评分参考序列后，盲用状态和候选字节保持完全一致。Python独立评分复现了所有原生正确、错误和非设计索引输出计数。','',
        '| 决策规则 | 模拟运行点 | 出现错误输出的点 | 错误输出次数 | 最大单点错误数 |',
        '|---|---:|---:|---:|---:|']
    for name,label in names.items():
        z=sim[sim.policy==name];lines.append(f'| {label} | {len(z)} | {int((z.false_accepted>0).sum())} | {int(z.false_accepted.sum())} | {int(z.false_accepted.max())} |')
    lines+=['','各条件共享设计或计数图，尺度之间还是嵌套的；上表的错误次数不是独立试验计数，也不能直接换算为未来样本的误差概率。参数扫描需要额外计算，不能继承单设置速度表中的优势倍数。','',
        'CRC数量目标有一个直接的确定性关系。设选中点接受A条，其中F条为错误或非设计输出，正确数为S=A−F；以真值挑选的网格最佳正确数为S*。最大化A保证A≥S*，所以 **0≤S*−S≤F**。'
        '零观察错误时，选中点必然达到相应单设置网格的最佳正确数；这一关系不保证F=0。所有模拟和真实点均验证此关系。'
        '因此，CRC计数可以提供不依赖参考真值的目标，但其可靠性仍取决于错误接受受到控制，不能把这条关系当作新的准确率保证。','',
        '## 2 规模揭示小文库中看不见的错误','',
        '以下为相同规则分段间距45 nt、目标保留碱基深度100×、内部回收q=0.05、独立替换率1%的配对规模实验。每个数字为三个独立随机家族的均值；±为家族间样本标准差。每个家族的三个规模共享前缀设计及对应随机流。','',
        '| k | 设计数 | 规则 | 正确恢复率 % | 错误输出数 |',
        '|---:|---:|---|---:|---:|']
    for k in [21,31]:
        for n in [1024,8192,65536]:
            for name in ['default','global_crc_r','union_consensus','two_support_consensus']:
                z=scale.loc[(k,n,name)]
                lines.append(f'| {k} | {n:,} | {names[name]} | {z[("correct_pct","mean")]:.3f} ± {z[("correct_pct","std")]:.3f} | {z[("false_accepted","mean")]:.2f} ± {z[("false_accepted","std")]:.2f} |')
    traces=tr[tr.n==65536]
    longdonor=int((traces.longest_other_match>=40).sum())
    lines+=['',f'对全部原生参数点中不同的错误候选进行了事后追踪：65,536设计规模下有 **{len(traces)}** 个不同的“运行图、索引、错误序列”组合，其中 **{longdonor}** 个含有至少40 nt与其他设计完全相同的片段。'
        '这些是图路径混合的直接线索，而非把错误简单解释为输出中的少数碱基错误。追踪未参与候选筛选。'
        '同一错误序列在多个参数点出现，不会被当成多个独立错误机制。完整索引、候选、核苷酸距离、供体区间和真路径可接受性均保留在源表中。','',
        '一致候选并集要求所有成功设置的完整序列相同，但某个错误候选可以在多个相关设置中重复出现，而其他设置仅失败或歧义；因此“一致”不等于“正确”。'
        '至少两个r支持也不是两次独立证实。这个失败例必须保留，不能只展示并集增加了多少恢复。','',
        '## 3 相同平均错误率可以有不同的恢复边界','',
        '匹配间距80 nt、q=1、目标30×/100×，比较独立替换、平均错误串长6 nt的平稳成簇替换，以及每条设计的固定模板位置替换。'
        '后者在该设计的所有分子中保持同样的错误；它是持续性损伤/系统性错误的极端对照，不是已拟合的化学衰变模型。三者预期每碱基替换率均为3%，实现误差另行记录。','',
        '| k | 目标深度 | 错误结构 | 默认正确恢复率 % |','|---:|---:|---|---:|']
    for (k,D,channel),row in error.iterrows():lines.append(f'| {k} | {D}× | {channel} | {row["mean"]:.3f} ± {row["std"]:.3f} |')
    lines+=['','增加深度可提高独立错误下某个正确窗口被观察到的机会，却不能保证恢复所有拷贝中都被固定改变的窗口。'
        '错误率必须连同相关结构一起描述；(1−e)^k不能直接用于持续性错误或插入缺失。','',
        '## 4 从经验相图推进到可检验的模型边界','',
        '独立Bernoulli断裂时，位于j的k窗口被完整保留的概率为 '
        '**w_j=(1−1/L)^(k−1)[q+(1−q)(A_j+B_j−A_jB_j)]**，其中A_j=(1−1/L)^j，B_j=(1−1/L)^(200−j−k)。'
        '规则分段则枚举L种相位得到精确平均。端点片段总保留，内部片段以q回收。','',
        '令每条设计的物理分子数服从Poisson(μ)，μ按期望保留碱基预算设定；独立替换下，正确来源窗口计数的边际分布为Poisson(λ_j)，λ_j=μw_j(1−e)^k。'
        '定义f_j=P(C_j≤c)，则完整输出路径支持的概率满足 **max(0,1−Σf_j)≤P(全部窗口支持)≤1−max f_j**。'
        '此处没有假设相邻窗口独立；直接相乘仅作为忽略依赖的对照。已知种子和未报告终端look-ahead仍需单独处理。','',
        f'在 **{v["theoretical_conditions"]}** 个条件中，**{v["theory_intersections"]}** 个条件的经验支持率，其Bonferroni校正Clopper–Pearson区间与理论概率范围相交。'
        '每个条件合并三个独立文库的3,072条设计，族错误水平0.05。这里的“相交”表示未检出不相容，不是证明模型正确，也不是宣称区间很窄。','',
        f'另做条件化丰度扩展：显式使用模拟生成的每条设计丰度，而不把丰度混合视为同一Poisson分布。'
        f'采用适用于独立但非同分布支持指示量的Hoeffding同时区间，**{v["conditional_abundance_intersections"]}/{v["conditional_abundance_conditions"]}** 个条件与条件化边界相交。'
        '这部分是测量后推导的无拟合参数扩展，单独标注，不能冒称全部预注册。','',
        '正确来源支持与原生图计数不是同一个对象。后者还可能包含跨设计碰撞或错误生成的匹配词。所有检查中，原生计数不低于正确来源计数（按原生饱和上限截断后比较）。'
        '模型边界只约束指定生成模型下的正确来源输出覆盖概率，不直接约束所有解码器、实际化学样本或完整文件恢复。','',
        '## 5 真实176天文库复核','',
        '使用原始双端FASTQ，不合并、不裁剪，k=31、c=1/3/5、八个r。285,583个查询索引全部进入盲选目标，包括不属于210,000条设计的索引；未用设计成员列表预先过滤目标。'
        '与首轮重合的16个参数点，所有非时间诊断值一致。这个文库已在第一轮观察过，应称外部行为复核，不是新的独立实验文库。','',
        '| 规则 | 正确输出 | 错误及非设计输出 | 选定 c / r |','|---|---:|---:|---|']
    for name,label in names.items():
        row=real[real.policy==name].iloc[0]
        setting='逐索引合并' if pd.isna(row.selected_r) else f'{int(row.selected_c)} / {row.selected_r:g}'
        lines.append(f'| {label} | {int(row.correct):,} | {int(row.false_accepted)} | {setting} |')
    lines+=['','## 6 稿件应如何完善','',
        '将“接近覆盖上限”限定为参考辅助、实现相关的精确分解：输出覆盖→任意终端支持→种子支持→原生整数窗口允许→唯一正确输出。'
        '原先连续r_needed和固定参考look-ahead的路径统计只能作为保守/描述性指标，不能冒充实际整数规则的最小r。旧图中的原测量数值保留，相关文字和图注必须明确口径。','',
        '本轮可加入模型、错误结构和规模/选参三组图。主文保留真实样本速度证据，同时明确多参数扫描没有做新的等资源速度比较；不能把候选选参包装为已经内置于生产DBGPS2的自动功能。'
        '实验局限包括随机载荷、有限规模、非标定的损伤模型、同一真实文库的重复利用，以及没有新的字节级完整文件恢复。','',
        '下一阶段最有价值的证据是独立实验文库的前瞻验证、测序深度与文件恢复概率的联系、强内码下的错误接受控制，以及在共同总时间预算下验证自适应选择。'
        '当前结果加强了可推广机制的证据，但不能替代这些实验或保证期刊层级。','',
        '图件均有PNG、PDF和SVG：Fig14_effective_coverage_model、Fig15_error_structure、Fig16_blind_selection_scale。'
        '完整锁定设计、候选字节、正确来源支持、原生计数和逐索引错误追踪均在交付数据包中。']
    lines += ['', '## 7 利用首次跌破阈值收紧边界','',
        '检查边际概率相乘的偏差后，进一步推导了不要求全路径Markov性的边界。相邻两个正确来源窗口的计数可分成共享Poisson成分和两个私有Poisson成分；共享均值来自完整(k+1)窗口的保留概率和错误存活率。由此精确计算 d_j=P(C_j>c, C_{j+1}≤c)。任意失败路径要么首窗口失败，要么至少有一次向下穿越阈值，因此 **P(全部支持)≥max(0,1−f_0−Σd_j)**。这仍是概率下界，不是独立窗口乘积或拟合曲线。','',
        f'新边界在原验证及全新4,096设计文库中核验：**{cv["compatible_conditions"]}/{cv["homogeneous_conditions"]}** 个条件的同时区间与边界相容；其中全新确认条件为 **{cv["fresh_compatible_conditions"]}/{cv["fresh_conditions"]}**。族水平0.05，采用精确Clopper–Pearson及Bonferroni调整。条件概率范围宽度中位数从 **{100*cv["median_original_interval_width"]:.2f}** 缩至 **{100*cv["median_crossing_interval_width"]:.2f}** 个百分点；某些失败过渡区域仍然很宽。','',
        '具体反例：Bernoulli间距80 nt、深度10×、q=0.05、e=1%、k=21、c=1时，原验证的完整来源支持为85.64%；边际独立乘积仅6.92%，窗口并集下界为0，而首次跌破阈值下界为77.94%，上界96.90%。这个结果把“窗口高度相关”的描述推进为可计算且可检验的范围。推导发生在初始结果之后，新确认种子在测量前另行锁定；两者不能混写为同一次预注册。','',
        '图14使用收紧后的边界，并保留独立乘积的失败例。模型只针对指定Poisson分子和独立替换过程，不把它外推成所有真实DNA样本的普适上限。']
    (R/'REPORT.md').write_text('\n'.join(lines)+'\n')
    key=dict(validation=v,crossing_validation=cv,total_new_graphs=366,total_new_operating_points=8784,scale65536={str(k):{name:{'correct_mean':float(scale.loc[(k,65536,name)][('correct_pct','mean')]),
            'correct_sd':float(scale.loc[(k,65536,name)][('correct_pct','std')]),
            'false_total':int(s[(s.k==k)&(s.n==65536)&(s.policy==name)].false_accepted.sum())} for name in names} for k in [21,31]},
            real={name:dict(correct=int(real[real.policy==name].iloc[0].correct),false_accepted=int(real[real.policy==name].iloc[0].false_accepted)) for name in names},
            donor_traces_65536=len(traces),long_donor_matches_65536=longdonor)
    (R/'key_findings.json').write_text(json.dumps(key,indent=2));print('DEEP_REPORT_COMPLETE',flush=True)


if __name__=='__main__':main()
