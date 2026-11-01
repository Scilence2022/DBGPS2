#!/usr/bin/env python3
"""Build bilingual, evidence-linked reports from frozen experiment records."""
import csv
import hashlib
import html
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'input/antkowiak2020'
RESULTS = DATA / 'results'
OUT = ROOT / 'output/antkowiak2020_reports'
PDF = ROOT / 'output/pdf/antkowiak2020'


def read(name): return json.loads((RESULTS / name).read_text())
def pct(n, d): return f'{n / d:.4%}'
def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


class Report:
    def __init__(self, number, language):
        self.number = number; self.zh = language == 'zh-CN'; self.language = language; self.lines = []
    def p(self, zh, en): self.lines.extend([zh if self.zh else en, ''])
    def h(self, zh, en): self.p('## ' + zh, '## ' + en)
    def table(self, zh, en, rows):
        headers = zh if self.zh else en
        self.lines.extend(['| ' + ' | '.join(map(str, headers)) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'])
        self.lines.extend('| ' + ' | '.join(str(v).replace('|', '/') for v in row) + ' |' for row in rows)
        self.lines.append('')
    def code(self, value): self.lines.extend(['```sh', value, '```', ''])


def common(report, summary, grid):
    r = report; f = r.number; j = summary[f'file{f}']; total = j['graph']['truth']
    raw_name = {1: 'I16_S2_R1_001.fastq.gz', 2: 'I18_S3_R1_001.fastq.gz', 3: 'File3.fq.gz'}[f]
    raw = json.loads((DATA / 'raw' / (raw_name + '.verified.json')).read_text())
    r.p(f'# File {f} 完整解码实验报告', f'# File {f}: complete decoding experiment report')
    r.p('版本：2026-10-18；包括本次会话全部已记录阶段。数值由本地 JSON/TSV 实验记录生成，英文版与中文版使用同一份数据。',
        'Version: 2026-10-18. Includes all recorded stages in this task. Numerical tables are generated from local JSON/TSV records; both languages use the same data.')
    r.h('1. 结论与实验范围', '1. Outcome and scope')
    if f < 3:
        r.p(f'File {f} 已从测序 reads 生成的实际片段输出完成外码解码。最佳严格图输出正确 {j["graph"]["correct"]:,} 条、错误 {j["graph"]["errors"]} 条，最终 ZIP 通过逐字节和内部 CRC 校验。未使用设计真值补片。',
            f'File {f} was decoded from actual strand decisions generated from sequencing reads. The best strict graph output contains {j["graph"]["correct"]:,} correct strands and {j["graph"]["errors"]} wrong strands. The recovered ZIP passes byte equality and internal CRC validation. No design strands were injected.')
    else:
        r.p('File 3 尚未从现有测序结果独立恢复完整文件。参考辅助分块后，前 11 块有满足外码界的实际输出；第 12 块仍未达标。使用设计真值补片的两种模拟均完整解码成功，必须与测序恢复区分。',
            'No complete File 3 has been recovered from the measured sequencing outputs alone. After reference-assisted read partitioning, actual outputs for blocks 1-11 satisfy the outer-code bound; block 12 does not. Two explicitly design-supplemented simulations do recover the complete file, and are reported separately.')
    r.p('本文的 DBGPS2 指实验适配器：复用本仓库精确 canonical k-mer 计数、图查询和动态覆盖比规则，增加论文所需的中部索引双向搜索及 GF(64) 校验。不是原版 CLI 对该数据格式的直接支持。没有把参考 payload 用作图搜索种子；File 3 后续 reads 归组步骤单独标明参考辅助。',
        'DBGPS2 here means an experimental adapter using this repository\'s exact canonical k-mer counter, graph queries and dynamic coverage-ratio rule, with an internal-index bidirectional traversal and GF(64) validation. It is not direct support by the stock CLI. Reference payloads do not seed graph traversal. Reference assistance in the later File 3 read-routing step is disclosed separately.')
    r.h('2. 数据来源、完整性与预处理', '2. Inputs, integrity and preprocessing')
    r.table(['项目', '值'], ['Item', 'Value'], [
        ['FASTQ', raw_name], ['Figshare file ID', raw['file_id']], ['Compressed bytes', raw['bytes']],
        ['MD5', raw['md5']], ['SHA-256', raw['sha256']], ['Raw reads', j['reads']['total_reads']],
        ['Retained reads', j['reads']['filtered_reads']], ['Design strands', total], ['Outer blocks', total // 16383]])
    r.p('全部原始文件通过公布的 MD5。基础图实验与同输入 MSA 使用全部长度 55-70 nt 的 reads，并将大于 60 nt 的序列截至前 60 nt；55-59 nt reads 保留实际长度。不是抽样分析。预处理元数据中的 sample_reads 仅对应额外生成的试运行文件，正式统计不使用它。',
        'The complete compressed input passed its published MD5. Baseline graph and matched MSA runs use every read of original length 55-70 nt; reads longer than 60 nt are capped at the first 60 nt, and shorter retained reads keep their length. These are not subsampled experiments. The sample_reads metadata describes an additional pilot file that is not used for the reported full-data results.')
    if f == 1:
        r.p('File 1 的物理设计序列和原 ZIP 来自作者仓库。后续公开设计副本与原设计一致。完整参考重新编码对照精确重现全部 16,383 条物理码字。',
            'Physical File 1 designs and the original ZIP come from the authors\' repository. The later public design copy agrees with the primary design. Re-encoding reproduces all 16,383 physical strands exactly.')
    else:
        r.p('File 2/3 设计来自后续公开数据仓库，其存储碱基方向和前端索引布局不等于物理合成格式。prepare_references.py 反转碱基顺序（不是反向互补），保留 42 nt payload，按原掩码重排中部索引并重新编码内码校验。该转换在 File 1 的全部码字上验证，转换后的本文件参考能生成 CRC 完整的 ZIP。',
            'Later File 2/3 design deposits do not use the physical strand orientation and index layout. prepare_references.py reverses base order, not reverse complements, preserves the 42-nt payload, restores the scrambled internal index and regenerates inner parity. This transformation is validated on every primary File 1 strand, and the transformed reference for this file yields a ZIP with valid CRCs.')
        validation = read('reference_raw_validation.json')
        r.p(f'在前 1,000,000 条原始 reads 的参考格式对照中，找到 {validation[f"File{f}_first_1000000_exact_reads"]} 条完整精确匹配。这里的百万条仅是参考格式验证，正式重建仍使用完整保留读集。完整文件比较基准由完整正确设计解码得到，而非从缺失的测序结果假定得到。',
            f'A reference-format check on the first 1,000,000 raw reads finds {validation[f"File{f}_first_1000000_exact_reads"]} full exact read matches. That subset is only a layout validation; reconstruction uses the complete retained input. The byte-comparison baseline is decoded from the complete validated design, not assumed to have been recovered from incomplete sequencing output.')
    r.h('3. 编码结构与外码恢复条件', '3. Encoding and outer-code criterion')
    r.p('物理布局为 42 nt payload + 12 nt scrambled index + 6 nt inner parity，共 60 nt；A/C/G/T 对应 0/1/2/3。索引为 24 bit，异或掩码 1071068。内码 RS(20,18) 使用 GF(2^6)、多项式 91，可修正一个 6-bit 符号。它不是通用 indel 纠错器；图搜索先生成全长候选，再进行符号纠错。',
        'The physical layout is 42-nt payload + 12-nt scrambled index + 6-nt inner parity, totaling 60 nt. A/C/G/T map to 0/1/2/3. The 24-bit index is XOR-scrambled with 1071068. Inner RS(20,18) over GF(2^6), polynomial 91, corrects one 6-bit symbol. It is not a general indel decoder: traversal produces full-length candidates before symbol correction.')
    r.p('外码为每块 6 路独立 RS(16383,10977)/GF(2^14)，冗余 5,406 个符号。每条片段携带 6 个 14-bit 外码符号。每一路须满足 2E_j + S <= 5406；S 为缺失片段数，E_j 为该路实际错误符号数。无错误时至少需要 10,977 个不同 index。完整片段错误数给出保守上界，后续实验按六路真实错误分别核算。各大块独立编码，没有可跨块补偿的额外外码。',
        'Each block contains six independent RS(16383,10977) codewords over GF(2^14), with 5,406 parity symbols. A strand carries one 14-bit symbol in each lane. Every lane must satisfy 2E_j + S <= 5406, where S is missing strands and E_j is wrong symbols in that lane. With no errors, 10,977 distinct indices suffice. Wrong-strand counts give a conservative bound; later results use actual lane errors. Blocks are independent, with no additional cross-block parity.')
    r.h('4. 图搜索规则、统计口径与对照', '4. Traversal, metric definitions and controls')
    r.p('只枚举图中实际观察到、包含有效 12 nt 索引的种子；在中部索引端点位置 54 锚定，向左延伸到起点，再补右端 6 nt 校验。每代最多 1050 条活跃路径，溢出时整个源 index 搜索视为擦除。全长严格候选须通过内码校验且唯一。k-mer 计数饱和于 16383，规范化时合并反向互补。',
        'Seeds are observed graph k-mers containing a valid 12-nt index. Traversal anchors at internal-index end position 54, extends left to the strand start, then right through the final 6 nt. At most 1,050 paths remain active per generation; overflow erases the entire source-index search. Strict output requires a unique full-length inner-code-valid strand. Canonical counts combine reverse complements and saturate at 16,383.')
    r.p('参数 c 的条件是 count > c，而不是 >= c。对上一步覆盖 p，相邻覆盖 q 必须位于 [max(c+1, floor(p/r)), min(16383, floor(p*r))]。所有正式网格 c>=1；大数据的两遍 singleton 过滤不改变可用 count>=2 键。只把 c 改成 0 并不能恢复单次支持，因为种子枚举也要求 count>1。',
        'Cutoff c means count > c, not >= c. For previous coverage p, next coverage q must be within [max(c+1, floor(p/r)), min(16383, floor(p*r))]. All measured grids use c>=1; the two-pass singleton filter preserves eligible counts >=2. Changing c alone to zero would not retain singletons because seed enumeration also requires count>1.')
    r.p('四种量分别报告：(1) 原始/预处理 reads 中直接观察的正确片段；(2) 候选池中至少出现一次完整正确 60 nt 的 index 数；(3) 每 index 最多一条实际输出中的正确、错误和擦除；(4) 完整文件字节及 ZIP 校验。跨参数候选并集要用真值标记成功，不能当作一个已实现的无参考选择器。',
        'Four outcomes are kept separate: (1) exact strands directly observed in input reads; (2) indices with at least one exact full 60-nt candidate; (3) correct, wrong and erased strands after an actual one-sequence-per-index decision; and (4) complete-file bytes and ZIP validation. Cross-parameter candidate unions require truth labels and are not an implemented reference-free selector.')
    if f == 1:
        r.p('原版 CLI 在实测和完美参考诊断输入上均输出 0 条；其前端索引、固定引物和 CRC/字节 RS 语义与本格式不同。随后排查随机数分布历史版本和外码频域/系统表示，确认使用历史 libc++ 均匀采样以及频域消息解码。正对照：完美参考 k=21 恢复全部 16383 条并生成原 ZIP。最佳真实图输出在不提供参考 FASTA 时逐字节一致。',
            'Stock-CLI diagnostic runs emit zero strands on both real and perfect-reference inputs because their prefix-index, primer and CRC/byte-RS semantics do not match this format. Compatibility diagnosis established the historical libc++ uniform sampler and spectral-message outer decoding. A perfect-reference k=21 positive control recovers all 16,383 strands and the original ZIP. The best real graph output is byte-identical when no evaluation FASTA is supplied.')
    if f == 2:
        r.p('早期 File 2 诊断曾直接使用后续设计的存储方向，也测试过前端索引解释。这些结果保存在 orientation_diagnostic，未纳入最终 35 组图网格；它们不是新的失败恢复率。最终物理格式对照、原始 reads 匹配和参考 ZIP 校验共同支持正确转换。',
            'Early File 2 diagnostics used the later deposited base order directly and also explored a prefix-index interpretation. Those artifacts remain under orientation_diagnostic and are excluded from the final 35-setting grid; they are not final recovery rates. Physical-layout checks, exact raw-read matches and the reference ZIP jointly validate the corrected conversion.')
    r.h('5. 未拆分 reads 的 MSA 与图解码尝试', '5. Pooled-read MSA and graph attempts')
    r.p('同输入 MSA 使用 prefix=16、最小簇 5、MUSCLE 3.8.31、gap 权重 0.4；较大簇取 5 次最多 15 条的随机子集，固定种子 20261017。它是 prefix 分簇基线，未完整复现论文 LSH。运行时间受并发与 worker 数影响，仅作为本机记录。',
        'Matched-input MSA uses a 16-base prefix, minimum cluster size 5, MUSCLE 3.8.31 and gap weight 0.4. Large clusters produce five subsets of at most 15 reads with fixed seed 20261017. This is a prefix-clustering baseline, not a full reproduction of the paper\'s LSH procedure. Recorded timings reflect concurrent workloads and differing workers, not controlled speed comparisons.')
    msa = j['msa']; inner = j['msa_inner']; graph = j['graph']
    r.table(['方法', '正确 index', '错误输出', '备注'], ['Method', 'Correct indices', 'Wrong outputs', 'Interpretation'], [
        ['Direct exact input', msa['raw_exact_recovered'], '-', pct(msa['raw_exact_recovered'], total)],
        ['MSA candidate coverage', msa['recovered'], '-', pct(msa['recovered'], total)],
        ['MSA strict decisions', inner['strict']['correct'], inner['strict']['errors'], 'GF(64) check'],
        ['MSA corrected decisions', inner['corrected']['correct'], inner['corrected']['errors'], 'One-symbol correction'],
        ['Graph strict decisions', graph['correct'], graph['errors'], f"k{graph['k']}/c{graph['cutoff']}/r{graph['ratio']}"]])
    r.p(f'MSA 产生 {msa["clusters"]:,} 个簇、{msa["candidates_generated"]:,} 条候选，其中 {msa["candidates_60nt"]:,} 条长度为 60 nt。完整处理时间 {msa["total_seconds"]:.3f} 秒；此数不是跨算法受控速度结论。内码纠错会提高正确输出，也可能产生误纠正，须同时看错误数。',
        f'MSA processes {msa["clusters"]:,} clusters and generates {msa["candidates_generated"]:,} candidates, of which {msa["candidates_60nt"]:,} are 60 nt. Recorded total processing time is {msa["total_seconds"]:.3f} s, not a controlled cross-algorithm benchmark. Inner correction can improve correct decisions while also introducing miscorrections, so wrong outputs are reported alongside gains.')
    if f == 1:
        n = j['notebook']; ni = j['notebook_inner']
        r.p('另按作者 notebook 的预处理设置运行：prefix=14、最小簇 4、gap 权重 0.5，先对长度 55-70 nt 的完整 reads 比对，投票后再截至 60 nt。这与同输入基线存在预处理差异，单独列示。',
            'An additional run follows the released notebook: prefix 14, minimum cluster 4, gap weight 0.5, alignment of complete 55-70-nt reads before capping the voted candidate at 60 nt. Its preprocessing differs from the matched baseline and is reported separately.')
        r.table(['notebook 方案', '正确', '错误'], ['Notebook variant', 'Correct', 'Wrong'], [
            ['Candidate coverage', n['recovered'], '-'], ['Strict decisions', ni['strict']['correct'], ni['strict']['errors']],
            ['Corrected decisions', ni['corrected']['correct'], ni['corrected']['errors']]])
        r.p(f'该 notebook 严格输出也独立完成 File 1 外码解码，ZIP 字节一致且 CRC 通过；实际用时 {j["notebook_decode"]["seconds"]:.3f} 秒。',
            f'The notebook strict output also independently decodes File 1 to byte-identical ZIP content with valid CRCs; recorded outer-decode time is {j["notebook_decode"]["seconds"]:.3f} s.')
    r.p(f'最终基础图参数表共 {len(grid)} 组；每个 k 的最佳实际严格输出如下。完整逐参数结果及所有保存的运行命令、退出码和时间见随附 evidence/File{f}/graph_parameters.csv 和 run_ledger.json。',
        f'The final pooled strict grid contains {len(grid)} settings. The best strict output at each k is listed below. All individual settings and saved commands, exit codes and timings are included in evidence/File{f}/graph_parameters.csv and run_ledger.json.')
    rows = []
    for k in sorted({a['k'] for a in grid}):
        subset = [a for a in grid if a['k'] == k]; b = max(subset, key=lambda x: (x['correct'], -x['errors']))
        rows.append([k, len(subset), b['cutoff'], b['ratio'], b['correct'], b['errors'], b['overflow']])
    r.table(['k', '参数数', '最佳c', '最佳r', '正确', '错误', '溢出index'], ['k', 'Settings', 'Best c', 'Best r', 'Correct', 'Wrong', 'Overflow indices'], rows)
    r.p('小 k 提高局部覆盖，同时增加跨链共享和含错图分支；较大 k 则损失低支持路径。因此不能只根据真值路径覆盖判断实际恢复能力。最佳参数由同一数据真值事后筛选，属于探索性实验，没有独立留出集。',
        'Smaller k increases local support but also cross-strand sharing and noisy branches; larger k loses low-support paths. True-path support alone therefore does not establish actual recoverability. Best parameters are selected retrospectively against this dataset\'s truth, without a held-out validation set.')


def finish_small(r, j):
    f = r.number; g = j['graph']; d = j['decode']
    r.h('6. 分块拆分是否执行及逐块外码核验', '6. Read-partition status and per-block verification')
    if f == 1:
        r.p('File 1 只有一个外码块，因此不存在跨块 reads 拆分前后对比。本文件没有运行参考比对分块。完整 79 组图尝试均基于未按块拆分的保留 reads。不能把 File 3 的分块改善套用到 File 1。',
            'File 1 has one outer-code block, so there is no cross-block read-partition comparison. Reference-alignment read routing was not run for this file. All 79 graph settings use the retained pooled input; File 3 partitioning gains are not transferred to File 1.')
    else:
        r.p('File 2 有 3 个外码块。本任务没有对 File 2 执行按设计比对拆分 reads 的实验，故没有“拆分后”的测量值；现有混池结果已使 3/3 块满足外码界并恢复完整 ZIP。下表按输出 index 统计各块，这不等于先拆分 reads 再解码。',
            'File 2 has three outer-code blocks. Design-alignment read partitioning was not performed for File 2, so no post-partition measurement is available. The pooled result already satisfies the bound in all three blocks and recovers the full ZIP. The table groups decoded indices by block; it does not imply that reads were partitioned before graph decoding.')
    r.table(['块(1起)', '正确', '错误', '擦除', '保守2E+S', '达标'], ['Block (1-based)', 'Correct', 'Wrong', 'Erased', 'Conservative 2E+S', 'Pass'], [
        [b['block'] + 1, 16383 - b['erasures'] - b['errors'], b['errors'], b['erasures'], b['budget'], b['within_rs_bound']] for b in g['per_block']])
    r.p(f'全部 {len(g["per_block"])} 块达标。扫描中共有 {j["configurations_with_all_blocks_in_bound"]} 组参数的全部块满足保守界；这不是说每组都运行过最终外码，最终文件记录以实际 decode JSON 为准。',
        f'All {len(g["per_block"])} blocks pass. In the scan, {j["configurations_with_all_blocks_in_bound"]} settings satisfy the conservative bound in every block. This does not mean that each was fully decoded; actual file-decoding attempts are identified by the decode records.')
    r.h('7. 完整文件与实现核验', '7. Complete-file and implementation validation')
    r.table(['验证项', '结果'], ['Validation', 'Result'], [
        ['Outer decoder exit', d['exit_code']], ['Decoded bytes', d['bytes']], ['SHA-256', d['sha256']],
        ['Byte equality', d['byte_identical']], ['ZIP CRC error', d['zip_error']], ['Outer seconds', f"{d['seconds']:.6f}"]])
    r.p('rs_driver_fast 使用作者频域外码算法、历史随机化，并采用等价的 GF(2^14) 对数表加速。232072 项有限域代数检查以及 File 1 原始/加速外码的字节交叉验证保留在 evidence。文件成功以实际 ZIP 内容验证，不只依据片段恢复率。',
        'rs_driver_fast uses the authors\' spectral outer-code algorithm and historical randomization with an equivalent GF(2^14) log-table acceleration. Evidence includes 232,072 finite-field checks and original-versus-accelerated File 1 byte comparison. File success is established from the decoded ZIP, not fragment recall alone.')
    r.p('完整文件：artifacts/File' + str(f) + '_sequencing_recovered.zip。File 1 对照作者原 ZIP；File 2 对照完整物理参考解出的 CRC 有效 ZIP。',
        f'Full output: artifacts/File{f}_sequencing_recovered.zip. File 1 is compared to the authors\' original ZIP; File 2 is compared to the CRC-valid ZIP decoded from the complete physical reference.')


def file3(r, summary):
    ratio = read('file3_ratio_summary.json'); inner = read('inner_rs/summary.json'); union = read('inner_rs/parameter_union/summary.json')
    split = read('block_decode/summary.json'); small = read('block12_small_k/summary.json'); sim = read('file3_supplement_simulation/summary.json')
    r.h('6. 未拆分混池：小 k 与覆盖比扩展', '6. Before partitioning: smaller k and ratio expansion')
    r.p('原始固定 r=2.55 的 35 组配置最高为 k19/c1：正确 135575、错误 111。随后增加 144 组，合计 179 组，包含 k=18/20 和更小 k、r=1.2 到 10 的组合。最高正确严格输出成为 k19/c1/r3.5：137662 条正确、143 条错误；达标 8/12 块。按达标块数最优的是 k20/c1/r3.5，正确 136445、错误 61，达标 9/12。',
        'The original 35 settings at r=2.55 peak at k19/c1: 135,575 correct and 111 wrong strands. A further 144 settings expand the grid to 179, including k18/k20, smaller k, and ratios from 1.2 to 10. The best strict correct count is k19/c1/r3.5: 137,662 correct, 143 wrong, with 8/12 blocks in bound. Maximizing passing blocks instead selects k20/c1/r3.5: 136,445 correct, 61 wrong, and 9/12 blocks in bound.')
    r.table(['块', '正确', '错误', '擦除', '保守预算', '达标'], ['Block', 'Correct', 'Wrong', 'Erased', 'Conservative budget', 'Pass'], [
        [b['block']+1, 16383-b['erasures']-b['errors'], b['errors'], b['erasures'], b['budget'], b['within_rs_bound']] for b in summary['file3']['graph']['per_block']])
    r.p('第 12 块在这一最佳严格输出中只有 6193 条正确，22 条错误，缺失 10168，保守预算 10212。k=13/15 虽然提高了真值 k-mer 覆盖，但大量分支导致溢出。覆盖上界不能直接当作图输出。',
        'Block 12 has only 6,193 correct and 22 wrong strands in that best strict result, with 10,168 erasures and a conservative budget of 10,212. Although k13/k15 increase true-k-mer support, extensive branching causes overflow. Graph support is not equivalent to emitted reconstruction.')
    r.h('7. 候选级内码 RS 与全部参数并集', '7. Candidate-level inner RS and the complete parameter union')
    r.p('先对每个未溢出完整图路径做 RS 单符号纠错，而不是对已通过内码的严格输出再次纠错。解码后 index 可以重分配；越界 index 拒绝。11 种实际决策包括 strict、corrected_unique、corrected_vote1/2/3/5/10/20、pooled_weight1/4/16；票数是图路径数而非 reads 数。存在零修正候选时 priority 规则优先唯一零修正码字；其余按唯一性或票数决定，并列擦除。pooled 规则联合零/一修正票数加权选唯一最大。',
        'One-symbol RS correction is applied to every complete, non-overflow graph path, not reapplied to strict outputs that already pass parity. Repaired indices can be reassigned; out-of-domain indices are rejected. Eleven actual rules are evaluated: strict, corrected_unique, corrected_vote1/2/3/5/10/20, and pooled_weight1/4/16. Votes count graph paths, not reads. Priority rules prefer a unique zero-correction word; other decisions use uniqueness or vote thresholds, with ties erased. Pooled rules choose a unique maximum using combined weighted zero/one-correction votes.')
    r.p('第一阶段利用完整路径数量上界证明 119 组不可能超过已达到的单配置正确数，仅重跑 60 组、660 个决策输出。此剪枝只适用于单配置最大值，不适用于跨参数并集。因此第二阶段补跑全部 119 组，最终覆盖原始 179 组；所有严格输出与原版本逐字节核对。',
        'Initially a complete-path-count upper bound excludes 119 settings from improving the attained single-setting maximum, leaving 60 reruns and 660 decision outputs. This pruning is valid for a single-setting maximum, not for a union. All 119 omitted settings are therefore later run, completing the original 179. Strict outputs are checked byte-for-byte against the earlier version.')
    r.table(['统计量', '结果'], ['Metric', 'Result'], [
        ['Best single candidate coverage', 140797], ['Best actual unique correct', 138674], ['Wrong at that actual output', 17144],
        ['179-setting correct candidate union', union['candidate_union_correct']], ['Zero-correction union', union['zero_candidate_union_correct']],
        ['RS required across all settings', union['requires_rs_across_all_parameters']], ['All-policy output correct union', union['all_policy_output_union_correct']],
        ['Block12 candidate union', 6927]])
    r.p('单配置候选上限 k19/c1/r3.5 有 3361623 条完整路径：139101 条零修正有效路径、1000579 条一符号修正路径、2161363 条不可修正路径和 60580 条越界路径。正确候选为 140797，其中 2496 个 index 仅通过修正进入该参数候选池。最佳实际正确数来自 k19/c1/r3/pooled_weight16：138674 正确但 17144 错误，只有 5/12 块达标；完整性不能按正确数单独判断。',
        'At k19/c1/r3.5, 3,361,623 complete paths comprise 139,101 zero-correction valid paths, 1,000,579 one-symbol-corrected paths, 2,161,363 uncorrectable paths and 60,580 out-of-domain paths. Correct candidate coverage is 140,797, including 2,496 indices requiring correction at this setting. The best actual correct count is k19/c1/r3/pooled_weight16: 138,674 correct but 17,144 wrong, passing only 5/12 blocks. Correct count alone is not a file-recovery criterion.')
    r.table(['块', '179组候选并集', '未覆盖', '仅擦除假设达标'], ['Block', '179-setting candidate union', 'Missing', 'Oracle erasure bound'], [
        [b['block'], b['correct'], b['missing'], b['oracle_erasure_bound_met']] for b in union['per_block']])
    r.p('并集按纠错后 index 和完整 60 nt 真值计一次成功，总计 151202/196596（76.9100%）。前 11 块的真值筛选并集可满足纯擦除条件，第 12 块仍缺 9456。并集只给出候选覆盖上界，没有实现一个无需真值的跨参数择优算法，也没有在这一阶段恢复完整 File 3。',
        'The union counts each repaired decoded index once if its full exact 60-nt truth occurs under any setting: 151,202/196,596 (76.9100%). Truth-selected unions for blocks 1-11 satisfy an erasure-only criterion; block 12 still misses 9,456 strands. This is candidate coverage, not a reference-free cross-parameter selection algorithm or complete File 3 recovery.')
    r.h('8. DBGPS-analyzer 对第 12 块的诊断', '8. DBGPS-analyzer diagnosis of block 12')
    r.p('使用用户指定工具源码独立编译，未修改工具源码。13 次完整计数包含处理后 reads 的 9 个 k、原始 reads 的 k17/19/20（实际长度）及原始 k19 截至 60 nt。另用 DBGPS-links 对全库和 12 个块逐 k 分析。目标始终为完整 60 nt、无引物裁剪。',
        'The user-specified analyzer is compiled from its recorded source without changes. Thirteen full count runs cover nine k values on processed reads, k17/k19/k20 on original full-length reads, and original reads capped at 60 nt at k19. DBGPS-links analyzes the entire design and each of the twelve blocks at every k. Targets remain full 60-nt strands with no primer trimming.')
    r.table(['k19输入', '第12块全k-mer计数>1', '覆盖中位数'], ['k19 input', 'Block12 all k-mers >1', 'Median strand mean coverage'], [
        ['Processed 55-70/cap60', 6277, 4.12], ['All raw/cap60', 7768, 5.64], ['All raw/full length', 7901, 5.71]])
    r.p('k19/c1/r3.5 的 9925 个 RS 后失败 index 中，7430 有零覆盖真值 k-mer，2429 最弱仅单次，28 被覆盖比拒绝，38 有合格真值路径但搜索溢出。全库 k19 噪声/目标不同键数量比 Kn=27.516543，不是碱基错误率。第 12 块 GC=50.85%，k19 内部设计交叉链接 3（其他块 1-7），粗粒度设计特征未显示明显异常。池中 canonical k-mer 覆盖不等于目标专属 read 深度。',
        'Among 9,925 block-12 failures after RS at k19/c1/r3.5, 7,430 have an absent true k-mer, 2,429 have only singleton minimum support, 28 fail the ratio rule, and 38 have an eligible true path but overflow. The full-library k19 distinct-noise/target-key ratio Kn=27.516543 is not a base error rate. Block12 has 50.85% GC and three within-block design cross-links at k19 (other blocks: 1-7), without an obvious coarse design anomaly. Pooled canonical coverage is not target-assigned read depth.')
    r.p('9 个 k 的最小覆盖与原计数器逐 index 一致；全部 179 组参数的整数方向覆盖规则加溢出状态，准确预测零修正正确候选标记。这是诊断验证，不把可达真值路径直接算作已恢复输出。',
        'Per-index minimum coverage agrees with the original counter at all nine k values. Directed integer-rounded feasibility plus recorded overflow states exactly predicts zero-correction correct-candidate flags across all 179 settings. This validates the diagnosis, not a claim that all feasible truth paths were emitted.')
    r.h('9. 按设计比对归组 reads 与同读集重混对照', '9. Design-alignment read routing and the same-read re-pooling control')
    r.p('minimap2 2.31-r1302 对全部 195619515 条原始 reads 与全设计比对。核心参数：-a -x sr -k11 -w2 -n2 -m15 -s20 -A2 -B3 -O4,12 -E2,1 -N5 -p0.8 --secondary=yes --eqx -t12。归组用到完整 payload，属于参考辅助；只路由原 reads，不修正碱基或填补缺口。SAM 反向主比对恢复原始方向及质量顺序。',
        'minimap2 2.31-r1302 aligns all 195,619,515 raw reads to the complete design with -a -x sr -k11 -w2 -n2 -m15 -s20 -A2 -B3 -O4,12 -E2,1 -N5 -p0.8 --secondary=yes --eqx -t12. Full payloads assist assignment, so this is reference-assisted. Routing changes no read bases or gaps. Reverse-strand SAM primaries are restored to original sequence orientation and quality order.')
    r.p('broad 要求匹配>=20、query 对齐长度>=30、payload 重叠>=20、AS>=25、匹配比例>=60%；所有得分达到最佳合格得分80%的已报告合格比对必须同块。strict 另要求主比对 MAPQ>=10、匹配>=24、query>=35、payload>=25、AS>=40、比例>=70%。最多报告5条次优比对，无法保证未报告替代归属不存在。',
        'Broad assignment requires >=20 matches, >=30 aligned query bases, >=20 payload-overlap bases, AS>=25 and >=60% match fraction. Every reported qualifying alignment scoring at least 80% of the best qualifying hit must name the same block. Strict additionally requires primary MAPQ>=10, >=24 matches, >=35 query bases, >=25 payload bases, AS>=40 and >=70% match fraction. At most five secondary hits are reported; unreported alternatives may still exist.')
    r.p('broad 共 44581669 条，跨块歧义 15770 条，未归组 151022076 条。broad 12 块、歧义和未归组构成完整互斥读集；strict/inclusive 是重叠变体。全量 read name/sequence/quality 的 SHA-256 多重集计数、和及异或一致，验证全部原始序列与质量守恒；不是仅小样本检查。',
        'Broad pools contain 44,581,669 reads, with 15,770 cross-block ambiguous and 151,022,076 unassigned reads. Twelve broad pools plus ambiguous and unassigned reads form a complete disjoint partition; strict/inclusive variants overlap them. Full read-name/sequence/quality SHA-256 multiset totals, sums and XORs agree with the source, validating complete record conservation rather than only a pilot subset.')
    groups = {(a['mode'], a['block'], a['subset']): a for a in split['groups']}
    rows = []
    counts = split['split_counts']['blocks']
    for b in range(1, 13):
        rows.append([b, counts[b-1]['broad_reads'], groups[('historical_processed', b, 'common10')]['candidate_union_correct'],
                     groups[('assigned_pooled', b, 'common10')]['candidate_union_correct'], groups[('broad', b, 'common10')]['candidate_union_correct']])
    r.table(['块', 'broad reads', '旧处理混池', '同读集重混', '拆分独立'], ['Block', 'Broad reads', 'Old processed pool', 'Same reads re-pooled', 'Separate pools'], rows)
    r.p('共同10组参数：k15/c1/r1.8,2.55；k17/18/19 各 c1/r2.55,3.5；k20/c1/r3.5,5。新图计数保留实际长度，ANTK_READ_LIMIT=200。与原版 limit=60 的11个策略逐字节一致。所有新实验保留全局 index 域196596和1050路径上限，不将单块重新编号。',
        'The matched ten settings are k15/c1/r1.8,2.55; k17/18/19/c1/r2.55,3.5; and k20/c1/r3.5,5. New counting keeps actual read lengths with ANTK_READ_LIMIT=200. At limit=60, all eleven policies reproduce the earlier output bytes. Global index domain 196596 and the 1050-path cap are retained; block indices are never renumbered locally.')
    r.p('同批已分配 reads 重混后与分块使用完全相同的序列/质量，另有全记录多重集核验。第12块共同网格从重混8734增至分块9861（+1127），这才隔离拆分本身的影响；旧处理混池6739还同时涉及输入筛选与长度变化。核心实验完整完成319组；另有5组全原始未筛选混池辅助结果，后续因k18一次搜索耗时700.778秒且196322个index溢出而停止。未完成配置不记作0，也不纳入共同10组对照。',
        'Re-pooling and splitting use exactly the same assigned sequences/qualities, independently checked by a full-record multiset audit. In block12, matched-grid coverage rises from 8,734 re-pooled to 9,861 separated (+1,127), isolating separation itself. The older processed-pool count of 6,739 also changes filtering and length handling. All 319 core settings complete. Five all-raw unfiltered mixed auxiliary settings complete; further work stops after one k18 search takes 700.778 s and overflows 196,322 indices. Incomplete settings are neither zero outcomes nor part of the matched ten.')
    r.p('第12块 broad 扫描原179组后并集10033，对比旧混池6927，新增3111、丢失5，净增3106。共同10组 strict、inclusive 并集分别9544和9860，均未优于broad的9861。该阶段最佳实际外码预算输出为k17/c1/r3.5/corrected_vote2，9366正确、43错误、6974擦除，最差预算7060。',
        'Across the original 179 settings, separated broad block12 reaches 10,033 versus 6,927 in the old mixed pool: 3,111 gained and five lost, net +3,106. Strict and inclusive matched-ten unions are 9,544 and 9,860, neither above broad\'s 9,861. At this stage the lowest-budget actual block12 output is k17/c1/r3.5/corrected_vote2: 9,366 correct, 43 wrong, 6,974 erased, worst budget 7,060.')
    r.p('以下各块的实际输出按已测六路最大预算最小选取；选择参数使用真值评估，但文件内候选决策本身不读真值。第12块改列后续小k阶段的更新结果。',
        'Actual per-block outputs below minimize the measured worst-lane budget. Parameter selection is truth-informed exploratory evaluation; decisions within each file do not read truth. Block12 is updated to the later small-k result.')
    chosen = split['chosen_per_block'][:11] + [small['best_outer_budget']]
    r.table(['块', 'k/c/r/策略', '正确', '错误', '擦除', '最大预算'], ['Block', 'k/c/r/policy', 'Correct', 'Wrong', 'Erased', 'Max budget'], [
        [a['block'], f"{a['k']}/{a['cutoff']}/{a['ratio']}/{a['policy']}", a['correct'], a['errors'], a['erasures'], a['max_lane_budget']] for a in chosen])
    r.h('10. 第12块更小 k 与 r 的174组扩展', '10. Block12: 174 additional small-k and ratio settings')
    r.p('同一批1681943条broad reads上新增174组，和原179组合计353组。初始129组补上k14/16并扩展k13/15/17，随后20组细调k16与k15的r、25组细调k14阈值。除k/c/r外，输入、路径上限、索引域及RS决策均保持一致。完整新网格见block12_small_k/plan.json及block12_small_k/parameters.tsv；包含大量效果较差的参数，未删除失败。',
        'On the same 1,681,943 broad reads, 174 additional settings bring the block12 total to 353. An initial 129 add k14/k16 and extend k13/k15/k17; twenty refine k16/k15 ratios, and twenty-five refine k14 cutoffs. Input, path limit, index domain and RS policies remain unchanged apart from k/c/r. block12_small_k/plan.json and block12_small_k/parameters.tsv preserve the full grid, including poor outcomes.')
    r.table(['k', '参数数', '该k候选并集', '最佳单参数候选', '相对旧179新增'], ['k', 'Settings', 'Union at k', 'Best single candidate', 'New vs prior 179'], [
        [a['k'], a['parameters'], a['combined_union'], a['best_candidate']['candidate_correct'], a['newly_recovered_vs_all_baseline']] for a in small['by_k'] if a['k'] <= 17])
    r.p('最佳单参数为k16/c1/r3.5，正确候选10237，其中10029来自零修正路径，208在此参数必须内码纠错；r3.75也达到10237。353组跨参数并集10659/16383（65.0613%），相对旧10033新增626；205个index在全部参数中都必须RS修正才成功。仍缺5724，仅正确片段选出的纯擦除集合距离10977要求还差318。',
        'The best single candidate count is k16/c1/r3.5: 10,237, including 10,029 with a zero-correction path and 208 requiring inner correction at that setting; r3.75 ties at 10,237. The 353-setting union is 10,659/16,383 (65.0613%), gaining 626 over 10,033. Across all settings, 205 indices require RS correction. There are still 5,724 missing indices; a truth-selected correct-only set is 318 short of the 10,977 erasure-only requirement.')
    r.table(['实际策略(k16/c1/r3.5)', '正确', '错误', '擦除', '最差2E+S'], ['Actual policy (k16/c1/r3.5)', 'Correct', 'Wrong', 'Erased', 'Worst 2E+S'], [
        ['pooled_weight16', 10077, 1540, 4766, 7826], ['strict', 9973, 15, 6395, 6425]])
    r.p('较保守strict的六路错误为[13,14,15,15,15,15]，预算[6421,6423,6425,6425,6425,6425]。k13/14过短造成较多碰撞和溢出；最优k16不能简化为越小越好。对最佳预算参数完全不传评分参考重跑，11种输出及status/decisions/inner_stats均逐字节一致。分块过程本身仍使用参考。',
        'Strict lane errors are [13,14,15,15,15,15], with budgets [6421,6423,6425,6425,6425,6425]. Very short k13/k14 create more collisions and overflow; k16 being optimal does not mean smaller is always better. Re-running the lowest-budget setting without an evaluation reference yields identical bytes for all eleven policies and status/decisions/inner_stats. Earlier read routing remains reference-assisted.')
    r.h('11. 补入设计真值后的完整解码模拟', '11. Full-file decoding after explicit design supplementation')
    r.p('两种情形均保留前11块实际输出，仅在第12块缺失index补设计正确序列，按缺失index升序，不替换已有片段。情形A保留全部15条错误，需补1019条，得到10992正确+15错误；情形B先借助真值选出10659个已成功候选index，再补318条，得到10977正确、0错误。情形B的候选挑选和后续注入都涉及真值。',
        'Both simulations retain the actual outputs for blocks 1-11 and inject correct design strands only at missing block12 indices, in ascending index order, never replacing existing strands. Case A retains all 15 wrong strands and adds 1,019, yielding 10,992 correct +15 wrong. Case B first truth-selects the 10,659 successful candidate indices and adds 318, yielding 10,977 correct and no wrong strands. Both candidate selection and injection in case B use truth.')
    r.table(['情形', '补入', '正确', '错误', '擦除', '最大预算', 'ZIP字节一致'], ['Case', 'Injected', 'Correct', 'Wrong', 'Erased', 'Max budget', 'ZIP exact'], [
        [a['name'], a['injected_correct'], a['after']['correct'], a['after']['wrong'], a['after']['erasures'], a['after']['max_lane_budget'], a['decode']['byte_identical_to_reference']] for a in sim['cases']])
    r.p('实际运行原外码解码器处理全部12块，两种模拟均得到1268026字节的ZIP，逐字节等于完整设计解码基准，并通过所有条目CRC。SHA-256为a7d53ade64bab3f60252e75b28e1741d85fdd7b5b393ffddc95fd141466d8cbb。补片清单和原有错误均保留。少补1条时预算5407，仅计算超界，未实际运行该超界完整解码。',
        'The original outer decoder actually processes all twelve blocks in both simulations, producing a 1,268,026-byte ZIP identical to the complete-design baseline and passing every entry CRC. SHA-256 is a7d53ade64bab3f60252e75b28e1741d85fdd7b5b393ffddc95fd141466d8cbb. Injection ledgers and original wrong strands are retained. One fewer injected strand gives budget 5,407; that out-of-bound full decode is not run, only its budget is calculated.')
    r.p('结论限于：假设补入上述数量正确片段，外码实现可完整恢复文件。不能把补片模拟计作纯测序恢复，也不能将完整设计解码基准计作新实验恢复结果。',
        'The conclusion is conditional: with these explicitly injected correct strands, the outer implementation recovers the full file. Supplemented simulations and the complete-design baseline are not sequencing-only recovery results.')


def appendices(r, grid):
    f = r.number
    r.h('附录A. 完整网格与运行记录', 'Appendix A. Exact grid and execution records')
    r.p('下表按k与r列出实际存在于最终基础指标中的c，不代表未列组合也测试过。完整每组正确/错误、溢出、歧义、每块预算，以及命令/时间在证据目录提供。相同参数的早期运行可能被后续重跑覆盖，运行记录保留历史；最终指标每个参数只计一次。',
        'The table lists measured c values for each k/r present in final pooled metrics. Unlisted combinations are not implied to have been tested. Evidence contains individual correctness, errors, overflow, ambiguity, block budgets, commands and timings. Earlier repeats may be superseded by later runs; ledgers retain history, while final metrics count each setting once.')
    rows = []
    for k in sorted({a['k'] for a in grid}):
        for ratio in sorted({a['ratio'] for a in grid if a['k'] == k}):
            cuts = sorted({a['cutoff'] for a in grid if a['k'] == k and a['ratio'] == ratio})
            rows.append([k, ratio, ', '.join(map(str, cuts))])
    r.table(['k', 'r', '实际c列表'], ['k', 'r', 'Measured c values'], rows)
    r.h('附录B. 复现入口与脚本职责', 'Appendix B. Reproduction entry points and script roles')
    r.p(f'独立代码包：File{f}_code.zip。解压后保持目录结构，从包根目录运行。包含当时DBGPS2核心源码、实验适配器、便携RS头文件快照、MUSCLE源压缩包、固定版本信息与样例测试数据。原始reads不随代码包分发；fetch阶段显式下载对应完整FASTQ并做MD5。',
        f'Standalone code archive: File{f}_code.zip. Preserve its layout and run from its root. It contains the recorded DBGPS2 core, adapters, portable RS-header snapshot, MUSCLE source archive, pinned-version metadata and small validation fixtures. Raw reads are excluded; the explicit fetch stage downloads the matching complete FASTQ and verifies MD5.')
    r.code('python3 tools/reproduce.py build\npython3 tools/reproduce.py smoke\npython3 tools/reproduce.py fetch\npython3 tools/reproduce.py prepare\npython3 tools/reproduce.py graph\npython3 tools/reproduce.py msa')
    if f < 3: r.code('python3 tools/reproduce.py outer')
    else:
        r.code('python3 tools/reproduce.py inner\npython3 tools/reproduce.py partition\npython3 tools/reproduce.py small-k\npython3 tools/reproduce.py supplement')
        r.p('File 3 的inner阶段先需要基础179组status记录；partition阶段需要minimap2和内码构建；small-k需要已导出的broad块池；supplement需要small-k报告及前11块实际输出。各阶段既有依赖见代码包README。并集报告需要NumPy，图表另需Matplotlib，建议Python3.12。完整重跑包含数亿reads，运行成本显著。',
            'File3 inner analysis requires the original 179 status files; partition needs minimap2 and the inner build; small-k requires exported broad block reads; supplementation needs the small-k report and actual first-eleven-block outputs. Package README documents dependencies. Union reports need NumPy, plots additionally need Matplotlib; Python 3.12 is recommended. Full replay includes hundreds of millions of reads and substantial runtime.')
    roles = [
        ['download.py / fetch_sources.py', '校验下载、版本固定' if r.zh else 'Verified downloads and source pins'],
        ['prepare_reads.c / prepare_references.py', '长度预处理、物理格式参考转换' if r.zh else 'Length preprocessing and reference transcoding'],
        ['graph_adapter.c / run_graphs.py', '严格图搜索与参数网格' if r.zh else 'Strict traversal and parameter grid'],
        ['msa_baseline.py / inner_select.cpp', 'MSA与共同内码决策' if r.zh else 'MSA and common inner-code decisions'],
        ['rs_driver.cpp / fast_gf14.hpp / legacy_random.hpp', '外码解码、有限域加速和历史采样' if r.zh else 'Outer decode, field acceleration and legacy sampler']]
    if f == 3: roles += [
        ['graph_inner_rs.cpp / run_inner_graph.py', '候选单符号纠错及11种策略' if r.zh else 'Candidate correction and eleven policies'],
        ['run_block12_analyzer.py / report_block12_analyzer.py', '独立覆盖与图可达诊断' if r.zh else 'Independent coverage and graph diagnostics'],
        ['split_reads_by_alignment.py / split_aligned_blocks.cpp', '全量比对及不改序列的分块' if r.zh else 'Full alignment and sequence-preserving routing'],
        ['run_block_decode.py / report_block_decode.py', '同读集重混与分块对照' if r.zh else 'Same-read re-pooling and block comparison'],
        ['run_block12_small_k.py / report_block12_small_k.py', '174组扩展与353组并集' if r.zh else '174 added settings and 353-setting union'],
        ['simulate_file3_supplement.py', '显式真值补片和全文件核验' if r.zh else 'Explicit truth injection and full-file verification']]
    r.table(['脚本', '职责'], ['Script', 'Role'], roles)
    r.h('附录C. 限制、来源与证据定位', 'Appendix C. Limitations, sources and evidence map')
    r.p('所有最优参数均为本数据集探索性选择；没有独立留出或多随机种子置信区间。MSA不是完整LSH复现。覆盖/候选并集、实际唯一输出和完整文件成功不可互换。输入筛选、实际长度、分块及额外真值信息有分别注明。记录中的未尝试、中止和诊断性失败不伪装成已完成测量。',
        'All best settings are exploratory choices on this dataset, with no held-out validation or multi-seed confidence intervals. MSA is not a complete LSH reproduction. Coverage, candidate unions, actual decisions and file recovery are distinct outcomes. Filtering, length handling, partitioning and added truth information are disclosed separately. Unattempted, interrupted and diagnostic failures are not represented as completed measurements.')
    r.p(f'主要本地证据：evidence/File{f}/ 下的原始 summary.json、graph_metrics.json、graph_parameters.csv、run_ledger.json、输入verified记录及相关decode验证；MANIFEST.sha256登记全部交付文件摘要。每份报告及证据与原工作区路径的对应见EVIDENCE_MAP.json。',
        f'Primary local evidence is under evidence/File{f}/: original summary.json, graph_metrics.json, graph_parameters.csv, run_ledger.json, verified-input records and decode validation. MANIFEST.sha256 records delivered-file digests; EVIDENCE_MAP.json links each evidence item to its original workspace path.')
    r.p('[论文：DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3)。[作者仓库](https://github.com/MLI-lab/noisy_dna_data_storage)，固定41422be448b3c1da9041c3a2967d2754d2ad8b1a。[RS仓库](https://github.com/reinhardh/dna_rs_coding)，checkout固定455e1a5182bf64281f54d434a4ed5ca456c2373e，所用历史头文件8e7e3b5。[后续设计仓库](https://github.com/fml-ethz/dt4dds-challenges_notebooks)，固定861d3b95b52d1b121e9d85040f20c8ecf0cda2d4。',
        '[Article: DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3). [Primary repository](https://github.com/MLI-lab/noisy_dna_data_storage), pinned to 41422be448b3c1da9041c3a2967d2754d2ad8b1a. [RS repository](https://github.com/reinhardh/dna_rs_coding), checkout 455e1a5182bf64281f54d434a4ed5ca456c2373e, historical headers from 8e7e3b5. [Later designs](https://github.com/fml-ethz/dt4dds-challenges_notebooks), pinned to 861d3b95b52d1b121e9d85040f20c8ecf0cda2d4.')


def copy_evidence(f, evidence_map):
    folder = OUT / f'evidence/File{f}'; folder.mkdir(parents=True, exist_ok=True)
    common = ['summary.json', f'file{f}_graph_metrics.json', f'file{f}_reads.json', f'file{f}_msa.json', f'file{f}_msa_inner.json',
              f'file{f}_best_decode.json', 'source_pins.json', 'reference_transcoding.json', 'reference_outer_validation.json',
              'reference_raw_validation.json', 'reference_manifest.json', 'fast_gf_validation.txt', 'legacy_random_validation.txt']
    if f == 1: common += ['file1_notebook.json', 'file1_notebook_inner.json', 'file1_notebook_decode.json', 'reference_blind_check.json',
                          'positive_control.json', 'stock_real.log', 'stock_perfect_crc.log']
    if f == 2: common += ['orientation_diagnostic/README.txt']
    if f == 3:
        common += ['file3_before_ratio_sweep.json', 'file3_ratio_summary.json', 'file3_ratio_validation.json',
                   'inner_rs/summary.json', 'inner_rs/validation.json', 'inner_rs/file3_candidate_path_bounds.json',
                   'inner_rs/comparison.csv', 'inner_rs/parameter_union/summary.json', 'inner_rs/parameter_union/parameters.tsv',
                   'inner_rs/parameter_union/per_index.tsv', 'inner_rs/parameter_union/candidate_success_matrix.tsv.gz',
                   'inner_rs/parameter_union/validation.json', 'analyzer_block12/summary.json', 'analyzer_block12/provenance.json',
                   'analyzer_block12/validation.json', 'analyzer_block12/all_blocks.tsv', 'analyzer_block12/best_parameter_loss_breakdown.tsv',
                   'block_decode/summary.json', 'block_decode/validation.json', 'block_decode/split_provenance.json',
                   'block_decode/full_fastq_audit.json', 'block_decode/assigned_pool_audit.json', 'block_decode/raw_ancillary_scope.json',
                   'block_decode/candidate_comparison.tsv', 'block_decode/decision_comparison.tsv', 'block_decode/reads/split_counts.json',
                   'block12_small_k/summary.json', 'block12_small_k/plan.json', 'block12_small_k/validation.json',
                   'block12_small_k/independent_audit.json', 'block12_small_k/parameters.tsv', 'block12_small_k/per_index.tsv',
                   'block12_small_k/candidate_success_matrix.tsv.gz', 'block12_small_k/decision_metrics.tsv',
                   'file3_supplement_simulation/summary.json', 'file3_supplement_simulation/validation.json']
        for case in ['actual_plus_1019', 'oracle_union_plus_318']:
            common += [f'file3_supplement_simulation/{case}/injected_design_strands.tsv', f'file3_supplement_simulation/{case}/outer_decode.log', f'file3_supplement_simulation/{case}/result.json']
    for name in common:
        source = RESULTS / name
        if not source.exists(): raise FileNotFoundError(source)
        target_name = 'graph_metrics.json' if name == f'file{f}_graph_metrics.json' else name
        target = folder / target_name; target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target)
        evidence_map[str(target.relative_to(OUT))] = dict(source=str(source), sha256=digest(source))
    raw = {1:'I16_S2_R1_001.fastq.gz', 2:'I18_S3_R1_001.fastq.gz', 3:'File3.fq.gz'}[f]
    verified = DATA / 'raw' / (raw + '.verified.json')
    shutil.copy2(verified, folder / 'input_verified.json')
    evidence_map[str((folder / 'input_verified.json').relative_to(OUT))] = dict(source=str(verified), sha256=digest(verified))
    records = []
    for source in sorted(RESULTS.glob(f'file{f}_*.run.json')):
        records.append(dict(source=str(source), record=json.loads(source.read_text())))
    if f == 3:
        for sub in ['inner_rs', 'analyzer_block12', 'block_decode/decodes', 'block12_small_k/decodes']:
            for source in sorted((RESULTS / sub).rglob('*.run.json')):
                records.append(dict(source=str(source), record=json.loads(source.read_text())))
    (folder / 'run_ledger.json').write_text(json.dumps(records, indent=2))
    grid = read(f'file{f}_graph_metrics.json')
    fields = ['k', 'cutoff', 'ratio', 'emitted', 'correct', 'errors', 'erasures', 'exact_recall', 'precision', 'overflow', 'ambiguous', 'graph_supported', 'per_block']
    with (folder / 'graph_parameters.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(grid)
    artifacts = OUT / 'artifacts'; artifacts.mkdir(exist_ok=True)
    if f < 3: shutil.copy2(RESULTS / f'file{f}_best.zip', artifacts / f'File{f}_sequencing_recovered.zip')
    else:
        for case in ['actual_plus_1019', 'oracle_union_plus_318']:
            shutil.copy2(RESULTS / f'file3_supplement_simulation/{case}/simulated_recovered.zip', artifacts / f'File3_SIMULATED_{case}.zip')
    return grid


def make_pdf(markdown, target, number, language):
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether
    font = 'ReportUnicode'
    if font not in pdfmetrics.getRegisteredFontNames():
        pdfmetrics.registerFont(TTFont(font, '/System/Library/Fonts/Supplemental/Arial Unicode.ttf'))
        pdfmetrics.registerFontFamily(font, normal=font, bold=font, italic=font, boldItalic=font)
    width = A4[0] - 36 * mm
    normal = ParagraphStyle('Normal', fontName=font, fontSize=9.4, leading=14, spaceAfter=6, wordWrap='CJK', splitLongWords=True)
    title = ParagraphStyle('Title', parent=normal, fontSize=21, leading=28, textColor=colors.HexColor('#153b57'), spaceAfter=15)
    heading = ParagraphStyle('Heading', parent=normal, fontSize=13, leading=19, textColor=colors.HexColor('#153b57'), spaceBefore=11, spaceAfter=7, keepWithNext=True)
    cell = ParagraphStyle('Cell', parent=normal, fontSize=8, leading=11, spaceAfter=0)
    code = ParagraphStyle('Code', parent=normal, fontName='Courier', fontSize=7.4, leading=10.5, backColor=colors.HexColor('#f1f4f6'), borderPadding=6)
    def inline(value):
        value = html.escape(value)
        value = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)', lambda m: '<link href="' + m[2] + '" color="#14658b">' + m[1] + '</link>', value)
        value = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', value)
        return value.replace('`', '')
    story = []; lines = markdown.read_text().splitlines(); i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip(): i += 1; continue
        if line.startswith('```'):
            body = []; i += 1
            while i < len(lines) and not lines[i].startswith('```'): body.append(lines[i]); i += 1
            story.append(KeepTogether([Paragraph('<br/>'.join(html.escape(x) for x in body), code), Spacer(1, 6)])); i += 1; continue
        if line.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].startswith('|'):
                cells = [x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch('[-:]+', x) for x in cells): rows.append(cells)
                i += 1
            count = len(rows[0]); weights = [min(27, max(7, max(len(row[c]) for row in rows))) for c in range(count)]
            if count == 2: weights = [18, 47]
            widths = [width * w / sum(weights) for w in weights]
            values = [[Paragraph(inline(value), cell) for value in row] for row in rows]
            table = Table(values, colWidths=widths, repeatRows=1, hAlign='LEFT')
            table.setStyle(TableStyle([('BACKGROUND', (0,0),(-1,0),colors.HexColor('#dce8f0')),('VALIGN',(0,0),(-1,-1),'TOP'),
                                       ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')]),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#8aabbf')),
                                       ('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5)]))
            story.extend([table, Spacer(1,8)]); continue
        if line.startswith('# '): story.append(Paragraph(inline(line[2:]), title))
        elif line.startswith('## '): story.append(Paragraph(inline(line[3:]), heading))
        else: story.append(Paragraph(inline(line), normal))
        i += 1
    def page(canvas, doc):
        canvas.saveState(); canvas.setFont(font, 8); canvas.setFillColor(colors.HexColor('#526b7b'))
        canvas.drawString(18 * mm, A4[1] - 12 * mm, f'ANTKOWIAK 2020 / DBGPS2 / FILE {number} / {language}')
        canvas.line(18*mm, A4[1]-14*mm, A4[0]-18*mm, A4[1]-14*mm)
        canvas.drawString(18*mm, 12*mm, 'Experimental record | 2026-10-18')
        canvas.drawRightString(A4[0]-18*mm, 12*mm, str(doc.page)); canvas.restoreState()
    doc = SimpleDocTemplate(str(target), pagesize=A4, rightMargin=18*mm, leftMargin=18*mm, topMargin=22*mm, bottomMargin=21*mm,
                            title=f'File {number} decoding report ({language})', author='DBGPS2 experimental analysis')
    doc.build(story, onFirstPage=page, onLaterPages=page)


def main():
    OUT.mkdir(parents=True, exist_ok=True); PDF.mkdir(parents=True, exist_ok=True)
    reports = OUT / 'reports'; reports.mkdir(exist_ok=True)
    summary = read('summary.json'); evidence_map = {}
    for f in [1, 2, 3]:
        grid = copy_evidence(f, evidence_map)
        for language in ['zh-CN', 'en']:
            r = Report(f, language); common(r, summary, grid)
            if f < 3: finish_small(r, summary[f'file{f}'])
            else: file3(r, summary)
            appendices(r, grid)
            target = reports / f'File{f}.{language}.md'; target.write_text('\n'.join(r.lines))
            pdf = PDF / f'File{f}.{language}.pdf'; make_pdf(target, pdf, f, language); shutil.copy2(pdf, reports / pdf.name)
            print('REPORT COMPLETE', target, pdf, flush=True)
    (OUT / 'EVIDENCE_MAP.json').write_text(json.dumps(evidence_map, indent=2))


if __name__ == '__main__': main()
