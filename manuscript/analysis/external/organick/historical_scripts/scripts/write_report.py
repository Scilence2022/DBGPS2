#!/usr/bin/env python3
"""Write the Chinese results report from machine-readable measured outputs."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT / "results/summary.json").read_text())
rows = {row["run"]: row for row in summary}
pre = json.loads((ROOT / "results/preprocessing.json").read_text())
provenance = json.loads((ROOT / "results/provenance.json").read_text())

table = ["| 原始覆盖 | 方法 | 精确恢复条数 | 精确恢复率 | 错误且地址已知 | 错误且地址未知 | 未恢复的已知地址 | 输出正确率 |",
         "|---|---|---:|---:|---:|---:|---:|---:|"]
filtered_table = ["| 原始覆盖 | k-mer ≥2、ratio=2.55：过滤前错误输出 | 过滤后错误输出 | 过滤后精确恢复率 | 过滤后输出正确率 |",
                  "|---|---:|---:|---:|---:|"]
quality = ["| 样本 | FASTQ reads | 原始覆盖 | 发现前缀锚点的 reads | 可截取完整 A/C/G/T 编码区的 reads | 独立种子数 |",
           "|---|---:|---:|---:|---:|---:|"]
for sample in ["5x_raw", "10x_raw", "full"]:
    stats = pre["samples"][sample]
    coverage = f'{stats["mean_raw_coverage"]:.4f}×'
    quality.append(f'| {sample} | {stats["raw_reads"]:,} | {coverage} | {stats["anchored_reads"]:,} | {stats["complete_anchored_reads"]:,} | {stats["distinct_seeds"]:,} |')
    for cfg, label in [("c0_r2.55", "DBGPS2，k-mer ≥1，ratio=2.55"), ("c1_r2.55", "DBGPS2，k-mer ≥2，ratio=2.55"),
                       ("modal", "观测整段的唯一众数基线")]:
        row = rows[sample+"."+cfg]
        table.append(f'| {coverage} | {label} | {row["exact"]:,} | {row["exact_rate"]:.2%} | {row["wrong_known_address"]:,} | {row.get("outputs_unknown_address",0):,} | {row["missing"]:,} | {row["output_precision"]:.2%} |')
    raw = rows[sample+".c1_r2.55"]
    filtered = rows[sample+".c1_r2.55.nohp"]
    filtered_table.append(f'| {coverage} | {raw["wrong_known_address"]+raw.get("outputs_unknown_address",0):,} | {filtered["wrong_known_address"]+filtered.get("outputs_unknown_address",0):,} | {filtered["exact_rate"]:.2%} | {filtered["output_precision"]:.4%} |')

timing = []
for sample in ["5x_raw", "10x_raw", "full"]:
    info = json.loads((ROOT / "results" / (sample+".timing.json")).read_text())
    if info.get("wall_seconds") is not None:
        timing.append(f'- {sample}：建图及两种主要参数搜索总耗时 {info["wall_seconds"]:.1f} 秒。')
    else:
        timing.append(f'- {sample}：建图 CPU 时间 {info["graph_cpu_seconds"]:.2f} 秒；计数 ≥1 / ≥2 两种主要搜索 CPU 时间分别为 {info["search_cpu_seconds"]["c0_r2.55"]:.2f} / {info["search_cpu_seconds"]["c1_r2.55"]:.2f} 秒。')

full = rows["full.c1_r2.55"]
five = rows["5x_raw.c1_r2.55"]
report = f'''# nbt.4079 编码分析与 DBGPS2 实测

完成日期：2026-10-17。论文：[Organick et al., Random access in large-scale DNA data storage](https://www.nature.com/articles/nbt.4079)，Nature Biotechnology 36, 242–248 (2018)。

## 结论

**DBGPS2 的图构建和路径搜索可以用于该数据集的 DNA 片段恢复，但当前命令行不能直接按论文格式解码。** 我完成了原始测序数据、设计序列的下载及校验，并用独立适配器复用了未修改的 DBGPS2 核心。

在全部 607,150 条设计片段上，采用 k-mer 计数至少为 2、ratio=2.55 的设置：全量 {pre["samples"]["full"]["mean_raw_coverage"]:.2f}× 原始覆盖精确恢复 **{full["exact"]:,} 条（{full["exact_rate"]:.2%}）**；约 5× 原始覆盖精确恢复 **{five["exact"]:,} 条（{five["exact_rate"]:.2%}）**。这是完整 110 nt 编码区的逐碱基精确比较，并非仅比较种子之后的后缀。

这些结果证明了片段恢复层面的可用性；**没有证明完成原文件无误解码，也没有证明优于论文的 BMA + 外层 RS 完整流程**。

## 下载了什么

作者公开仓库为 [uwmisl/data-nbt17](https://github.com/uwmisl/data-nbt17)，本次固定版本为 `1077c74a0e002be83ab832e69bfd560e94d4d202`。

- 原始 FASTQ：`data/id20.fastq.gz`，1,363,110,383 字节，实测 **15,719,903 条 reads，均为 177 nt**。
- 原始参考文件：`data/id20.refs.txt.gz`，18,284,141 字节，包含 **607,150 条完整 150 nt 参考序列**。
- 标准 FASTA 转换：`data/references.fa`。保留所有参考序列，未改动碱基。
- 作者说明：`sources/MANIFEST`；LFS 指针、下载 URL、SHA256 校验记录均保存在 `sources/`。

两个压缩文件的 SHA256 均与作者固定版本的 Git LFS 记录一致；gzip 校验和与完整 FASTQ 格式检查通过。公开的文件名虽然叫 `id20`，MANIFEST 和补充说明将其指认为论文的 **File 16**：原始数字文件 9,503,066 字节，607,150 条片段。MANIFEST 中“逐 read 参考 ID”的描述与实际 `.refs` 内容不同，后者实际是一条设计序列一行。未将它误当成逐 read 标签，也没有用真值给 reads 分组。

该公开仓库只提供这一个文件的数据，未提供全部 35 个文件的原始数据，也未提供论文中的 nanopore 数据。本次结论适用于这里实测的 Illumina File 16。

## 论文如何编码

根据正文第 3 页、图 2 和补充说明：

1. **随机化**：数字文件与伪随机序列做 XOR，削弱重复模式和片段间相似性，并用于减少引物与数据区的非特异配对。
2. **跨片段 RS 外码**：分块后形成矩阵，每格为 16 bit；每块最多约 10 行、55,000 个数据列。对每一行做长 RS 编码，增加冗余列，每一列成为一条 DNA 片段。RS 码长可达 65,536。File 16 采用相对数据量 15% 的逻辑冗余，相当于全部编码片段约 13.04% 是冗余。
3. **地址和 DNA 内码**：加入 block/column 地址。每 3 个四进制符号表示 6 bit，转换成 4 个三进制符号，再通过依赖前一个碱基的旋转码映射为 DNA，避免同聚物。该变换的标称信息率为 1.5 bit/nt，之后还要扣除地址、外码、引物开销。
4. **文件级随机访问**：两端加文件专用引物。File 16 实测布局为：

```text
AGTGCAACAAGTCAATCCGT | 110 nt 编码区 | AATTGAATGCTTGCTTGCCG
       20 nt                           20 nt
```

正文第 4 页描述的恢复过程为：相似性聚类 → 修改过的 BMA 处理替换/插入/删除、重建片段 → 逆旋转及内码解码 → 外层 RS 修正错误和丢失 → 逆随机化恢复数字文件。

## 为什么需要适配

当前 DBGPS2 有三处格式假设不适用：默认引物是另一条固定 18 nt 序列；`-M` 只支持上下文无关的固定块映射，不能直接表达旋转码；默认片段内 CRC16 和现有 `-E rs:NB` 均不是论文的跨片段、16 bit 符号 RS 外码。因此，选择 `-E rs:16B` 不能代替论文的纠错层。

不过，de Bruijn 图直接处理 DNA k-mer，并不要求先将 DNA 逆编码成文件字节。适配器 `scripts/dbgps2_adapter.c` 直接调用原始 `count_file` 与 `find_paths_protected`，关闭不适用的片段内校验，以实际 reads 中的引物后前缀作为种子。生产源文件未修改，版本为 `{provenance["git_commit"]}`。

具体设置为 k=31，内部 index=16、path=96：输出 112 nt 中，前 2 nt 是真实引物尾部，后 110 nt 是完整编码区，去掉前 2 nt 即得到结果。这里的“16 nt index”只是接口布局；前 14 nt 在参考集中互不重复是本数据的实测特征，不能据此断言论文的地址字段恰好长 14 nt。**参考序列只用于评分，未用于建图、发现种子或修正输出。**

对照验证：20 条无误设计序列输入后，三种搜索配置均恢复 20/20；再对每条片段构造 5 条各含一个不同位置替换错误的 read、不给出完整无误 read，k-mer ≥2 设置仍恢复 20/20。原生 CLI 使用其固定引物时，对对应的无误输入不能恢复目标，验证了格式适配的必要性。这些对照不计入真实数据结果。

## 真实数据实验

固定随机种子 `20261017`，对原始 FASTQ 做均匀、嵌套抽样，没有按 reference 均衡覆盖。所有选中的原始 reads 均进入图，包括含 N 或不完整的 reads；两种方向都扫描锚点以发现种子。不使用真值聚类、不使用真值种子、不对错误 read 先做参考纠正。

{chr(10).join(quality)}

`-c 0` 表示 k-mer 计数至少为 1，`-c 1` 表示至少为 2。众数基线只是相同 reads 中观测到的完整 110 nt 序列的唯一众数，平票不输出；**它不是作者的 BMA 实现**。另尝试过关闭 ratio 的附加诊断，但搜索明显变慢，在主要配置完成后停止；其不完整输出以 `.tsv.incomplete` 保存，**不计入任何恢复率**。默认复现脚本仅运行已完成的两种主要配置。

{chr(10).join(table)}

“精确恢复率”以全部 607,150 条设计片段为分母；“输出正确率”以实际返回结果为分母，包含地址未知的错误输出。“未恢复”不含已输出错误序列的已知地址。不能仅看输出条数判断恢复质量。

![片段恢复与错误输出比较]({ROOT}/results/recovery_comparison.png)

所有参考编码区都没有相邻相同碱基。因此还测试了返回结果的无同聚物过滤；这只是论文内码的必要条件，并未实现完整内码解码。该过滤不能挽回核心搜索因多路径歧义而拒绝的片段。

{chr(10).join(filtered_table)}

## 与论文的比较边界

补充表 2（补充 PDF 第 4 页）中，**File 16 的最小平均可用 reads 覆盖为 5.0×，片段覆盖标准差为 3.0**。正文报告的是完整流程经过外层 RS 后，数字文件与原文件零字节差异；它并没有报告“5× 时 100% 的 DNA 片段都逐碱基恢复”。

本实验约 5× 是原始 FASTQ records / 607,150；其中可截取完整 A/C/G/T 编码区的 reads 约为 {pre["samples"]["5x_raw"]["complete_anchored_reads"]/607150:.2f}×。这不是论文原始可用 reads 过滤器的复现。我们只做了单个固定随机种子的抽样，未测出最低可成功解码覆盖。**原始/可用覆盖、片段/文件指标以及内外码校验能力都不相同，不能直接得出覆盖需求更高或更低的结论。**

可将 `(未恢复 + 2×已知地址错误)/607150` 作为保守的片段级错误/丢失粗略指标，机器结果保存在 `summary.tsv`；但实际 RS 约束逐 block、逐行计算，不同 DNA 错误影响的符号行也不同。即使该粗略值低于 13.04%，也不足以证明能纠回文件。

要完成真正的端到端比较，还需要作者完整内码/地址解释、外层 RS 参数和分块信息、随机化种子及原始数字文件。正文 Code availability（第 9 页）明确说明自定义代码是 proprietary，公开替代说明主要针对聚类。当前仓库没有提供这些完整解码材料，因此本次没有声称完成原始 PDF 的字节级恢复，也没有声称跑出了作者 BMA 的同条件结果。

## 建议

DBGPS2 适合作为“从测序 reads 重建候选 DNA 片段”的替代模块。若要将其用于论文的完整恢复流程，应保留论文的外层 RS，并为路径候选增加正确的旋转码/三进制块约束和地址解析，而不是使用不对应的 CRC 或片段内 RS。尤其在低覆盖下，允许计数为 1 的路径与抑制错误/歧义之间需要权衡。

## 复现与产物

```sh
python3 experiments/nbt4079/scripts/download_data.py
python3 experiments/nbt4079/scripts/count_stream.py
python3 experiments/nbt4079/scripts/prepare_reads.py --total 15719903
python3 experiments/nbt4079/scripts/run_experiment.py
python3 experiments/nbt4079/scripts/export_fasta.py
uv run --with matplotlib python experiments/nbt4079/scripts/plot_results.py
python3 experiments/nbt4079/scripts/write_report.py
```

- `results/summary.tsv`、`results/summary.json`：全部配置及无同聚物过滤结果。
- `results/*.tsv`：每条实际输出及其 read 种子和种子支持度。
- `results/*.c1_r2.55.nohp.fa.gz`：可直接使用的预测 DNA 片段 FASTA（110 nt 编码区）。保留无同聚物过滤后仍存在的错误，未用参考真值净化或修补。
- `results/preprocessing.json`：实际样本数、锚点统计、随机种子。
- `results/provenance.json`：源码版本、源文件 SHA256、编译命令和参数。
- `results/*.log`、`*.timing.json`：日志与运行耗时。
- `METHODS.md`：完整英文方法与比较口径。

预处理耗时 {pre["elapsed_seconds"]:.1f} 秒，以下仅是本机实际耗时，不与论文缺失的运行时间作性能排名：

{chr(10).join(timing)}

没有创建 commit、PR 或推送任何文件。下载数据和临时文件由本实验目录的 `.gitignore` 排除。
'''
(ROOT / "REPORT.zh-CN.md").write_text(report)
print(ROOT / "REPORT.zh-CN.md")
