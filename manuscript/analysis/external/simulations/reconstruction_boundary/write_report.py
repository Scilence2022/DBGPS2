"""Write measured, bounded conclusions; do not modify the Word manuscript."""
from pathlib import Path
import json
import pandas as pd
R=Path(__file__).resolve().parent
def main():
 d=pd.read_csv(R/'all_simulation_results.csv');v=json.loads((R/'validation_summary.json').read_text());real=pd.read_csv(R/'HT176_exact_boundaries.csv')
 default=d[(d.c==1)&(d.r==2.55)]
 def cell(cohort,L,D,q,e,k):return default[(default.cohort==cohort)&(default.length==L)&(default.depth==D)&(default.recovery==q)&(default.error==e)&(default.k==k)]
 def mean(s,f):return s[f].mean()
 def fmt(s,f):return f'{s[f].mean():.2f} ± {s[f].std():.2f}'
 biased=cell('main',45,100,.05,.01,21);confirm=cell('independent_confirmation',45,100,.05,.01,21)
 regular=cell('geometry_matched_regular',25,100,1,.01,31);geo=cell('geometric_breaks',25,100,1,.01,31)
 r255=real[(real.c==1)&(real.r==2.55)].iloc[0];r3=real[(real.c==1)&(real.r==3)].iloc[0];r4=real[(real.c==1)&(real.r==4)].iloc[0]
 table='| r | 正确重建 | 精确规则允许 | 覆盖上限 | 比例剪枝损失 | 搜索/歧义损失 |\n|---|---:|---:|---:|---:|---:|\n'
 for _,z in real[real.c==1].iterrows():table+=f'| {"off (10⁶)" if z.r==0 else f"{z.r:g}"} | {int(z.correct):,} | {int(z.ceil_adm):,} | {int(z.ceil_cov):,} | {int(z.ratio_losses):,} | {int(z.search_gap):,} |\n'
 text=f'''# DBGPS2 重建边界：首轮模拟与真实数据探索

模拟应当成为机制探索的主线，真实文库承担外部核验。首轮已实际完成 **{v['synthetic_graphs']:,} 个模拟计数图、{v['synthetic_operating_points']:,} 个参数运行点，以及 16 个 HT176 原始双端数据运行点**。参数点不是独立重复：每个条件使用三个独立随机文库。主矩阵每库 128 条设计，独立确认每库 512 条，几何/混合错误/丰度压力测试每库 256 条。全部运行点的损失分解及 {v['native_search_comparisons']:,} 次原生搜索对照通过。详细模型、先后锁定记录和限制见 PROTOCOL.md。

## 1. 当前结果支持什么

### 覆盖上限可以很高，但固定 r 仍会剪掉真路径

在规则分段间距45 nt、目标保留碱基覆盖100×、内部片段回收率q=0.05、独立替换错误率1%、k=21时，默认c=1、r=2.55的正确恢复为 **{fmt(biased,'correct_pct')}%**，覆盖上限为 **{fmt(biased,'ceil_cov_pct')}%**，平均差距 **{mean(biased,'ceil_cov_pct')-mean(biased,'correct_pct'):.2f} 个百分点**。独立512条设计文库在相同条件下正确恢复 **{fmt(confirm,'correct_pct')}%**、覆盖上限 **{fmt(confirm,'ceil_cov_pct')}%**。这里的±为三个独立模拟文库之间的样本标准差。

这说明默认2.55不是普适最优参数。高平均覆盖不能保证每一步相邻真k-mer计数都落入该窗口；必须检查局部过渡的尾部，而不能只看文库平均覆盖或全局CV。

### 平均碎片长度不能单独决定边界

在相同256条设计、目标100×、q=1、替换错误率1%、k=31以及相同名义切割密度1/25下，规则分段与Bernoulli随机断裂的观察平均片段长度分别为 **{fmt(regular,'mean_fragment_length')} nt** 和 **{fmt(geo,'mean_fragment_length')} nt**。两种模型的覆盖上限分别为 **{fmt(regular,'ceil_cov_pct')}%** 与 **{fmt(geo,'ceil_cov_pct')}%**，正确恢复分别为 **{fmt(regular,'correct_pct')}%** 与 **{fmt(geo,'correct_pct')}%**。

规则分段最大长度25 nt，不能提供31-mer；随机断裂具有较长片段的尾部，即使平均长度相近仍能提供31-mer。因此，“平均长度小于k就不能重建”是错误概括。可推广的对象应是覆盖每个连续k窗口的概率及其最弱位置。

### r 存在受数据约束的可用窗口

HT176原始双端数据中，c=1时的覆盖上限为 **{int(r255.ceil_cov):,}/{int(r255.n_ref):,} = {r255.ceil_cov_pct:.4f}%**。默认r=2.55正确恢复 **{int(r255.correct):,}** 条，距上限 **{int(r255.ceil_cov-r255.correct)}** 条（**{r255.ceil_cov_pct-r255.correct_pct:.4f}个百分点**）。r=3时差距缩至 **{int(r3.ceil_cov-r3.correct)}** 条；r=4虽然允许更多真路径，但 **{int(r4.overflow_adm):,}** 条规则本来允许的真路径被搜索膨胀阻断。

{table}

上表另有1条固定终端支持损失，已单独记录在完整CSV中。它使用未合并的原始双端文件，与前一轮成功合并分子的速度比较分开。可重合的12个旧全量原始数据运行点，其各类输出计数全部复现。这里只证明该真实文库存在窗口，不把2—3外推成所有文库的推荐范围。

## 2. 正确的上限分解

对固定图、k、c和r，定义 Ucov 为148个输出决定步骤均有至少c+1覆盖的设计数；Uend再要求某个终端延伸存在；Useed再要求当前实现有可用种子支持；Uadm再要求真路径逐步通过实际整数规则；S为唯一且正确的输出。

**S ≤ Uadm ≤ Useed ≤ Uend ≤ Ucov**。

总损失可以精确拆为输出覆盖缺口、终端支持缺口、种子支持缺口、比例剪枝损失以及搜索/唯一性损失。所有模拟运行点及真实数据运行点中，Uadm−S都等于允许真路径中的路径溢出与歧义损失之和。这是对该实现的参考辅助诊断，不是所有解码器的信息论上限，也不是完整文件恢复上限。

本轮还修正了两个近似口径：原生整数截断使连续对称比例条件偏保守；终端look-ahead不在报告片段内，必须允许任何合格终端延伸。最初使用固定参考终端碱基的开发运行已隔离，所有主条件均重跑，未混入统计。

## 3. 如何陈述“接近上限”

应同时满足：Ucov本身足够高；种子/终端支持缺口很小；窗口不会排除太多真过渡；允许的错误分支不会触发路径上限或有效候选歧义。只满足一个条件不够。作为描述性口径，主矩阵中覆盖上限至少50%的 {v['default_nontrivial_main_contexts']} 个均值条件里，默认设置有 {v['default_nontrivial_contexts_within_2pp']} 个距上限不超过2个百分点。2个百分点是探索性报告尺度，不是新的普适常数。

默认c=1、r=2.55下，本轮模拟记录中错误输出总数为 {v['default_wrong_outputs']}；全参数网格中单个参数运行点最多有 {v['all_settings_max_wrong_per_library']} 个错误输出。因此不能把“零错误”推广到所有参数。主矩阵、确认文库、混合错误、丰度偏倚和断裂几何的完整结果均保留。

## 4. 候选机制与下一步

在独立替换、无跨链碰撞等简化条件下，某个窗口的有效清洁支持近似由局部采样量 × 跨越整个k窗口的概率 × (1−e)^k 决定。规则随机相位分段的内区窗口存活因子约为 max(1−(k−1)/L,0)，独立Bernoulli断裂则为(1−1/L)^(k−1)；完整未切割分子应单独处理。重叠窗口高度相关，不能将各窗口成功概率简单相乘后称为严格恢复率公式。回收偏倚和插入缺失还会改变这个近似。

优先下一步是：在更大独立文库、真实系统性错误模型和独立实验文库中核验边界；构建不使用参考序列的r选择特征；在相同文件恢复要求下验证是否减少测序量。当前Uadm依赖真实设计作诊断，不能用真值挑选最佳r后宣称盲测最优。512条设计的确认也不能替代百万设计规模的碰撞压力测试。

## 5. 论文如何使用这些结果

可以先将主线改为“有效k窗口覆盖建立信息可用性边界；局部覆盖约束在保留真过渡与控制分支之间建立可用窗口”。保留失败区域和默认参数的反例。现有连续r_needed应标为保守诊断指标，上限应说明种子、终端和阈值口径。此次没有自动把探索性结果插入Word，也没有改变此前速度表。

图件：Fig_boundary_phase_maps（恢复/上限相图）；Fig_boundary_ratio_window（比例窗口及损失分解）；Fig_boundary_fragment_geometry（匹配断裂几何与覆盖诊断）。均提供PNG、PDF、SVG和源表。相图为三个文库的均值，图中的200*表示完整分子；HT176为一个真实文库，未人为添加独立重复误差线。

文献背景：原始DBGPS已讨论真路径完整性，因此不能把“路径覆盖上限”本身作为全新贡献。[Song et al., 2022](https://www.nature.com/articles/s41467-022-33046-w)。回收流程偏倚的模型动机来自[Gimpel et al., 2024](https://pubs.rsc.org/en/content/articlehtml/2024/dd/d4dd00220b)，本模拟未拟合或复刻其数字孪生。
'''
 (R/'REPORT.md').write_text(text)
 (R/'key_findings.json').write_text(json.dumps(dict(bias_example=biased[['correct_pct','ceil_cov_pct','ratio_losses_pct']].mean().to_dict(),bias_confirmation=confirm[['correct_pct','ceil_cov_pct','ratio_losses_pct']].mean().to_dict(),regular_geometry=regular[['correct_pct','ceil_cov_pct','mean_fragment_length']].mean().to_dict(),bernoulli_geometry=geo[['correct_pct','ceil_cov_pct','mean_fragment_length']].mean().to_dict()),indent=2))
 print('REPORT_WRITTEN')
if __name__=='__main__':main()
