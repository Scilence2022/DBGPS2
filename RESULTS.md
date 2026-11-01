# DBGPS2 × DNA Fountain (Erlich & Zielinski, Science 2017) 片段恢复测试报告

Worktree: `DBGPS2-dnafountain`（分支 `dnafountain-test`）
数据：ENA ERR1816980（主池 MiSeq PE151；R1 完整 16,692,317 对，md5 与 ENA 一致；
R2 下载中断致中段损坏，清洗后 7,828,857 干净读对 = 全 run 的 ~47%）
论文代码：github.com/TeamErlich/dna-fountain；原始输入 md5（论文 README）：
`8651e90d3a013178b816b63fdbb94b9b`（2,146,816 B = 67,088 × 32 B）

## 1. 论文编码方法

- 2,146,816 B 文件 → 67,088 段 × 32 B；Luby 变换（robust soliton，c=0.025，δ=0.001）
- droplet = 4B 种子 + 32B 载荷 + 2B RS(38,36)（reedsolo：GF(256)、0x11d、α=2、fcr=0）
- 种子 = 32 位 Galois LFSR 顺序状态（1+x^25+x^26+x^30+x^32，初态 0b001010101）；
  设计消耗约前 577k 个状态、按生化筛选（同聚物 ≤3、GC 45–55%）取 72,000 条（通过率 ~12.5%）
- 字节→DNA：{00,01,10,11}→{A,C,G,T}，MSB 在前 —— **与 DBGPS2 内置映射完全一致**
- 寡核苷酸 = 26nt 接头 + 152nt 数据块 + 21nt 接头；ENA 的 FASTQ 已去接头：
  R1 = block[0:151]，rc(R2) = block[1:152]（实测 52.8% 读对精确满足）

## 2. 结构兼容性与适配

| 项 | 兼容性 | 处理 |
|---|---|---|
| 2bit/碱基映射 | ✅ 完全一致 | 无需处理 |
| index 16nt 在链首 | ✅ | `-i 16` |
| 载荷 128nt 定长 | ✅ | `-l 128` |
| RS(2B) 冗余 | ✅ 参数逐位一致（已交叉验证 C↔reedsolo） | `-E rs:2B` |
| reads 以种子开头（无左侧上下文） | ⚠️ | R1 预置 10nt 恒定前缀 + `-A CCTGCAGAGT`（k=21>i=16） |
| 32 位稀疏 LFSR 种子（非连续区间） | ⚠️ | 新增 `-I 种子列表文件` |
| 论文解码语义 max_hamming=0 | ⚠️ | `DBGPS_DETECT_ONLY=1`（RS 只检不纠） |
| DBGPS2 搜索的 3′ 单碱基前瞻 | ⚠️ | 去接头 reads 缺该碱基 → R2 读段前置 'A'（反补端补 'T'） |

## 3. 测试结果

### 3.1 全量数据（7.83M 干净读对，~47% run）

| 管线 | 有效 droplet | 说明 |
|---|---|---|
| 论文式（读对精确合并 + RS 严格过滤） | 67,033 | 模拟论文 README 流程 |
| **DBGPS2**（k21, -c 0, -r 6, detect-only） | **71,983** | 99.7% 观测种子（72,383） |
| 交集 | 66,968，**99.996% 序列完全一致** | 双向交叉验证 |
| DBGPS2 独有 | 5,015 | 无完美读对、靠图共识恢复 |

**端到端解码（仅用 DBGPS2 的 71,983 条 droplet）：
67,088/67,088 段全部解开，输出文件 md5 = 8651e90d3a013178b816b63fdbb94b9b，
与论文原始输入文件逐位一致 —— 零错误完整恢复 2.15 MB。**

对比论文：论文解码器从完整 33M reads 消耗 69,870 条 oligo 完成恢复；
DBGPS2 用 ~47% 的 reads 即达成（更高单条质量 + 恢复了论文式管线漏掉的 5,015 条）。

### 3.2 一 tile 当量（前 750,000 读对 = 1.5M reads）

| 管线 | 有效 droplet |
|---|---|
| 论文式 | 65,428 |
| DBGPS2（修复前瞻后） | 63,327（94.4% of K，低于 ~104% 解码阈值，无法独立全解码；并集 68.6k 亦不足） |

低覆盖（λ≈11）下 DBGPS2 的固定比值窗口机制弱于「等一条完美读对」策略；
λ≳50 后 DBGPS2 反超。

### 3.3 参数敏感性（全量数据）

- k=16：跨寡核苷酸 chimera（池内 ~7k 共享 16-mer）→ k≥21 消除
- -r 2.55（默认）：λ≈6 时 136 步存活率 ~0.02%（泊松噪声撞窗口）→ -r 6–8 最优
- detect-only vs 纠正模式：产出基本相同（71,983 vs 71,961）
- 覆盖率 -c：低覆盖下 -c 0 唯一可行；-c≥2 剪枝过狠

## 4. 关键发现（对 DBGPS2 项目的反馈价值）

1. **3′ 单碱基前瞻**：搜索循环 `while(p <= path_len)` 比名义路径多走一步，
   依赖 reads 含 3′ 侧翼碱基（原始 DBGPS 数据有引物）；去接头外部数据集会全军覆没
   （表现：搜索死亡）。建议增加「无前瞻」开关或文档说明。
2. **固定比值窗口 vs 泊松覆盖**：低覆盖下真路径必然被剪；-s 只保护链首交界。
3. **-I 稀疏索引列表**：外部数据集普遍使用非连续/随机 32 位索引，值得合入主线。
4. **-A 自定义引物**：同上（本分支已实现）。
5. 观察到一次同命令不同结果的非确定性（12,618 vs 23,809，复跑 4 次稳定 23,809），
   未复现，建议留意。
6. 朴素多数共识在嵌合亚群下不可靠（失配案例中 90% 是共识错、DBGPS2 对）；
   RS 码字校验是更强判据。

## 5. 复现命令

```bash
# 1) 种子提取（LFSR 校验）
python3 analysis/extract_seeds.py --r1 r1g.fa --out seeds_all.tsv
# 2) R1 预置前缀、R2 前置 'A'（补 3′ 前瞻碱基）
# 3) DBGPS2 片段恢复（detect-only 对齐论文 max_hamming=0）
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -I seedlist.txt -o out r1pg.fa r2pT.fa
# 4) 喷泉码解码（Python3 复刻含 Python2 random 语义）→ md5 校验
python3 analysis/decode_fountain.py out.strands 的序列列 --no_screen --out decoded.bin
```

## 6. 产物位置

- 重建链：`analysis/runs/full/dbgps_final_do.strands`（71,983 × 152nt）
- 解码文件：`analysis/runs/full/decoded_noscreen.bin`（md5 8651e90d…）
- 各规模运行目录：`analysis/runs/{onetile,tile750k,full}/`
