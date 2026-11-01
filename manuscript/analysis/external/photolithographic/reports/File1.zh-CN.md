# File 1 完整解码实验报告

版本：2026-10-18；包括本次会话全部已记录阶段。数值由本地 JSON/TSV 实验记录生成，英文版与中文版使用同一份数据。

## 1. 结论与实验范围

File 1 已从测序 reads 生成的实际片段输出完成外码解码。最佳严格图输出正确 16,125 条、错误 0 条，最终 ZIP 通过逐字节和内部 CRC 校验。未使用设计真值补片。

本文的 DBGPS2 指实验适配器：复用本仓库精确 canonical k-mer 计数、图查询和动态覆盖比规则，增加论文所需的中部索引双向搜索及 GF(64) 校验。不是原版 CLI 对该数据格式的直接支持。没有把参考 payload 用作图搜索种子；File 3 后续 reads 归组步骤单独标明参考辅助。

## 2. 数据来源、完整性与预处理

| 项目 | 值 |
| --- | --- |
| FASTQ | I16_S2_R1_001.fastq.gz |
| Figshare file ID | 22305708 |
| Compressed bytes | 994447128 |
| MD5 | 6d48dee09e7cbe3f11bc749c9a7839c6 |
| SHA-256 | 25b70f73fcb4582666541f5379a74dc86d38e5e4ec690d8fad002986b03dad48 |
| Raw reads | 29303855 |
| Retained reads | 14299688 |
| Design strands | 16383 |
| Outer blocks | 1 |

全部原始文件通过公布的 MD5。基础图实验与同输入 MSA 使用全部长度 55-70 nt 的 reads，并将大于 60 nt 的序列截至前 60 nt；55-59 nt reads 保留实际长度。不是抽样分析。预处理元数据中的 sample_reads 仅对应额外生成的试运行文件，正式统计不使用它。

File 1 的物理设计序列和原 ZIP 来自作者仓库。后续公开设计副本与原设计一致。完整参考重新编码对照精确重现全部 16,383 条物理码字。

## 3. 编码结构与外码恢复条件

物理布局为 42 nt payload + 12 nt scrambled index + 6 nt inner parity，共 60 nt；A/C/G/T 对应 0/1/2/3。索引为 24 bit，异或掩码 1071068。内码 RS(20,18) 使用 GF(2^6)、多项式 91，可修正一个 6-bit 符号。它不是通用 indel 纠错器；图搜索先生成全长候选，再进行符号纠错。

外码为每块 6 路独立 RS(16383,10977)/GF(2^14)，冗余 5,406 个符号。每条片段携带 6 个 14-bit 外码符号。每一路须满足 2E_j + S <= 5406；S 为缺失片段数，E_j 为该路实际错误符号数。无错误时至少需要 10,977 个不同 index。完整片段错误数给出保守上界，后续实验按六路真实错误分别核算。各大块独立编码，没有可跨块补偿的额外外码。

## 4. 图搜索规则、统计口径与对照

只枚举图中实际观察到、包含有效 12 nt 索引的种子；在中部索引端点位置 54 锚定，向左延伸到起点，再补右端 6 nt 校验。每代最多 1050 条活跃路径，溢出时整个源 index 搜索视为擦除。全长严格候选须通过内码校验且唯一。k-mer 计数饱和于 16383，规范化时合并反向互补。

参数 c 的条件是 count > c，而不是 >= c。对上一步覆盖 p，相邻覆盖 q 必须位于 [max(c+1, floor(p/r)), min(16383, floor(p*r))]。所有正式网格 c>=1；大数据的两遍 singleton 过滤不改变可用 count>=2 键。只把 c 改成 0 并不能恢复单次支持，因为种子枚举也要求 count>1。

四种量分别报告：(1) 原始/预处理 reads 中直接观察的正确片段；(2) 候选池中至少出现一次完整正确 60 nt 的 index 数；(3) 每 index 最多一条实际输出中的正确、错误和擦除；(4) 完整文件字节及 ZIP 校验。跨参数候选并集要用真值标记成功，不能当作一个已实现的无参考选择器。

原版 CLI 在实测和完美参考诊断输入上均输出 0 条；其前端索引、固定引物和 CRC/字节 RS 语义与本格式不同。随后排查随机数分布历史版本和外码频域/系统表示，确认使用历史 libc++ 均匀采样以及频域消息解码。正对照：完美参考 k=21 恢复全部 16383 条并生成原 ZIP。最佳真实图输出在不提供参考 FASTA 时逐字节一致。

## 5. 未拆分 reads 的 MSA 与图解码尝试

同输入 MSA 使用 prefix=16、最小簇 5、MUSCLE 3.8.31、gap 权重 0.4；较大簇取 5 次最多 15 条的随机子集，固定种子 20261017。它是 prefix 分簇基线，未完整复现论文 LSH。运行时间受并发与 worker 数影响，仅作为本机记录。

| 方法 | 正确 index | 错误输出 | 备注 |
| --- | --- | --- | --- |
| Direct exact input | 4439 | - | 27.0952% |
| MSA candidate coverage | 10779 | - | 65.7938% |
| MSA strict decisions | 10778 | 4 | GF(64) check |
| MSA corrected decisions | 13179 | 779 | One-symbol correction |
| Graph strict decisions | 16125 | 0 | k19/c1/r2.55 |

MSA 产生 215,090 个簇、350,470 条候选，其中 346,253 条长度为 60 nt。完整处理时间 281.007 秒；此数不是跨算法受控速度结论。内码纠错会提高正确输出，也可能产生误纠正，须同时看错误数。

另按作者 notebook 的预处理设置运行：prefix=14、最小簇 4、gap 权重 0.5，先对长度 55-70 nt 的完整 reads 比对，投票后再截至 60 nt。这与同输入基线存在预处理差异，单独列示。

| notebook 方案 | 正确 | 错误 |
| --- | --- | --- |
| Candidate coverage | 11938 | - |
| Strict decisions | 11934 | 3 |
| Corrected decisions | 13602 | 706 |

该 notebook 严格输出也独立完成 File 1 外码解码，ZIP 字节一致且 CRC 通过；实际用时 2.243 秒。

最终基础图参数表共 79 组；每个 k 的最佳实际严格输出如下。完整逐参数结果及所有保存的运行命令、退出码和时间见随附 evidence/File1/graph_parameters.csv 和 run_ledger.json。

| k | 参数数 | 最佳c | 最佳r | 正确 | 错误 | 溢出index |
| --- | --- | --- | --- | --- | --- | --- |
| 12 | 4 | 1 | 10.0 | 0 | 0 | 16383 |
| 13 | 9 | 199 | 2.55 | 65 | 1 | 7480 |
| 15 | 12 | 19 | 2.55 | 11091 | 1 | 3440 |
| 17 | 12 | 3 | 2.55 | 16047 | 0 | 24 |
| 19 | 10 | 1 | 2.55 | 16125 | 0 | 0 |
| 21 | 12 | 1 | 2.55 | 16006 | 0 | 0 |
| 25 | 10 | 1 | 2.55 | 15495 | 0 | 0 |
| 31 | 10 | 1 | 2.55 | 13983 | 0 | 0 |

小 k 提高局部覆盖，同时增加跨链共享和含错图分支；较大 k 则损失低支持路径。因此不能只根据真值路径覆盖判断实际恢复能力。最佳参数由同一数据真值事后筛选，属于探索性实验，没有独立留出集。

## 6. 分块拆分是否执行及逐块外码核验

File 1 只有一个外码块，因此不存在跨块 reads 拆分前后对比。本文件没有运行参考比对分块。完整 79 组图尝试均基于未按块拆分的保留 reads。不能把 File 3 的分块改善套用到 File 1。

| 块(1起) | 正确 | 错误 | 擦除 | 保守2E+S | 达标 |
| --- | --- | --- | --- | --- | --- |
| 1 | 16125 | 0 | 258 | 258 | True |

全部 1 块达标。扫描中共有 30 组参数的全部块满足保守界；这不是说每组都运行过最终外码，最终文件记录以实际 decode JSON 为准。

## 7. 完整文件与实现核验

| 验证项 | 结果 |
| --- | --- |
| Outer decoder exit | 0 |
| Decoded bytes | 99103 |
| SHA-256 | 05295bc085b243e16aafa18a63a7b5fad944037f68667a7ac60bff691ac15945 |
| Byte equality | True |
| ZIP CRC error | None |
| Outer seconds | 0.291609 |

rs_driver_fast 使用作者频域外码算法、历史随机化，并采用等价的 GF(2^14) 对数表加速。232072 项有限域代数检查以及 File 1 原始/加速外码的字节交叉验证保留在 evidence。文件成功以实际 ZIP 内容验证，不只依据片段恢复率。

完整文件：artifacts/File1_sequencing_recovered.zip。File 1 对照作者原 ZIP；File 2 对照完整物理参考解出的 CRC 有效 ZIP。

## 附录A. 完整网格与运行记录

下表按k与r列出实际存在于最终基础指标中的c，不代表未列组合也测试过。完整每组正确/错误、溢出、歧义、每块预算，以及命令/时间在证据目录提供。相同参数的早期运行可能被后续重跑覆盖，运行记录保留历史；最终指标每个参数只计一次。

| k | r | 实际c列表 |
| --- | --- | --- |
| 12 | 2.55 | 1, 3 |
| 12 | 10.0 | 1, 3 |
| 13 | 2.55 | 1, 3, 9, 19, 49, 99, 199 |
| 13 | 10.0 | 1, 3 |
| 15 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 15 | 10.0 | 1, 3 |
| 17 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 17 | 10.0 | 1, 3 |
| 19 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 21 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 21 | 10.0 | 1, 3 |
| 25 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 31 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |

## 附录B. 复现入口与脚本职责

独立代码包：File1_code.zip。解压后保持目录结构，从包根目录运行。包含当时DBGPS2核心源码、实验适配器、便携RS头文件快照、MUSCLE源压缩包、固定版本信息与样例测试数据。原始reads不随代码包分发；fetch阶段显式下载对应完整FASTQ并做MD5。

```sh
python3 tools/reproduce.py build
python3 tools/reproduce.py smoke
python3 tools/reproduce.py fetch
python3 tools/reproduce.py prepare
python3 tools/reproduce.py graph
python3 tools/reproduce.py msa
```

```sh
python3 tools/reproduce.py outer
```

| 脚本 | 职责 |
| --- | --- |
| download.py / fetch_sources.py | 校验下载、版本固定 |
| prepare_reads.c / prepare_references.py | 长度预处理、物理格式参考转换 |
| graph_adapter.c / run_graphs.py | 严格图搜索与参数网格 |
| msa_baseline.py / inner_select.cpp | MSA与共同内码决策 |
| rs_driver.cpp / fast_gf14.hpp / legacy_random.hpp | 外码解码、有限域加速和历史采样 |

## 附录C. 限制、来源与证据定位

所有最优参数均为本数据集探索性选择；没有独立留出或多随机种子置信区间。MSA不是完整LSH复现。覆盖/候选并集、实际唯一输出和完整文件成功不可互换。输入筛选、实际长度、分块及额外真值信息有分别注明。记录中的未尝试、中止和诊断性失败不伪装成已完成测量。

主要本地证据：evidence/File1/ 下的原始 summary.json、graph_metrics.json、graph_parameters.csv、run_ledger.json、输入verified记录及相关decode验证；MANIFEST.sha256登记全部交付文件摘要。每份报告及证据与原工作区路径的对应见EVIDENCE_MAP.json。

[论文：DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3)。[作者仓库](https://github.com/MLI-lab/noisy_dna_data_storage)，固定41422be448b3c1da9041c3a2967d2754d2ad8b1a。[RS仓库](https://github.com/reinhardh/dna_rs_coding)，checkout固定455e1a5182bf64281f54d434a4ed5ca456c2373e，所用历史头文件8e7e3b5。[后续设计仓库](https://github.com/fml-ethz/dt4dds-challenges_notebooks)，固定861d3b95b52d1b121e9d85040f20c8ecf0cda2d4。
