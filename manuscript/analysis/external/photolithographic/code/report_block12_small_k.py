#!/usr/bin/env python3
"""Score small-k outputs and index-wise unions against the frozen 179-run baseline."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import numpy as np
from report_block_decode import B, N, DATA, PAT, bits, nbits, score_outputs

OUT = DATA / 'results/block12_small_k'
BASE = DATA / 'results/block_decode'
START = 11 * B


def candidate_flags(stem):
    path = Path(stem + '.truth_candidates.tsv')
    with path.open('rb') as f:
        f.seek(max(0, path.stat().st_size - 256))
        assert f.read().splitlines()[-1].startswith(b'196595\t'), path
    with path.open() as f:
        next(f)
        rows = np.array([tuple(map(int, row.split())) for row in itertools.islice(f, START, N)], dtype=np.int64)
    assert rows.shape == (B, 4) and np.array_equal(rows[:, 0], np.arange(START, N))
    assert np.array_equal(rows[:, 3], (rows[:, 1] + rows[:, 2]) > 0)
    return bits(rows[:, 3]), bits(rows[:, 1] > 0), rows[:, 3].astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--partial', action='store_true')
    args = ap.parse_args()
    truth = [line for line in (DATA / 'sources/File3_synthesized.fasta').read_text().splitlines() if line and not line.startswith('>')]
    assert len(truth) == N
    plan = json.loads((OUT / 'plan.json').read_text())
    expected = {(r['k'], c, float(v)) for r in plan['plans'] for c, v in r['grid']}
    prior_scores = [r for r in json.loads((BASE / 'all_decision_metrics.json').read_text()) if r['mode'] == 'broad' and r['block'] == 12]
    scores = prior_scores[:]
    configs = []; matrices = []; baseline = new = zero = 0; by_k = {}; completed = set()
    sources = [('baseline', BASE / 'decodes/broad/block12')]
    sources.extend(('small_k', p) for p in sorted((OUT / 'decodes').glob('*/block12')))
    for experiment, folder in sources:
        for path in sorted(folder.glob('*.inner_stats.json')):
            match = PAT.fullmatch(path.name)
            k, c, r = int(match[1]), int(match[2]), float(match[3])
            stem = str(path)[:-len('.inner_stats.json')]
            if args.partial:
                try: cb, zb, flags = candidate_flags(stem)
                except (AssertionError, FileNotFoundError): continue
            else: cb, zb, flags = candidate_flags(stem)
            row = dict(experiment=experiment, k=k, cutoff=c, ratio=r, candidate_correct=nbits(cb), zero_candidate_correct=nbits(zb), candidate_path=stem + '.truth_candidates.tsv')
            configs.append(row); matrices.append(flags); zero |= zb
            by_k.setdefault(k, dict(baseline=0, small_k=0, combined=0, parameters=0))
            by_k[k][experiment] |= cb; by_k[k]['combined'] |= cb; by_k[k]['parameters'] += 1
            if experiment == 'baseline': baseline |= cb; continue
            completed.add((k, c, r)); new |= cb
            stats = json.loads(path.read_text())
            for policy, emitted in stats['emitted'].items():
                output = Path(stem + '.' + policy + '.strands')
                rows, correct = score_outputs(output, truth, [12])
                result = dict(experiment=experiment, mode='broad', k=k, cutoff=c, ratio=r, policy=policy, path=str(output), **rows[12])
                assert result['emitted'] + result['outside_target_block_outputs'] == emitted
                assert not (correct[12] & ~cb)
                scores.append(result)
            print('SCORED', k, c, r, row['candidate_correct'], flush=True)
    assert nbits(baseline) == 10033
    assert len([r for r in configs if r['experiment'] == 'baseline']) == 179
    if not args.partial: assert completed == expected, (len(completed), len(expected), expected - completed)
    combined = baseline | new
    best_candidate = max(configs, key=lambda r: r['candidate_correct'])
    best_unique = max(scores, key=lambda r: (r['correct'], -r['errors']))
    best_budget = min(scores, key=lambda r: (r['max_lane_budget'], -r['correct'], r['errors']))
    baseline_s = json.loads((BASE / 'summary.json').read_text())
    old_group = next(g for g in baseline_s['groups'] if g['mode'] == 'broad' and g['block'] == 12 and g['subset'] == 'all')
    summaries = []
    for k, group in sorted(by_k.items()):
        rows = [r for r in configs if r['k'] == k]; decisions = [r for r in scores if r['k'] == k]
        summaries.append(dict(k=k, parameters=group['parameters'], baseline_union=nbits(group['baseline']),
                              new_union=nbits(group['small_k']), combined_union=nbits(group['combined']),
                              newly_recovered_vs_all_baseline=nbits(group['small_k'] & ~baseline),
                              best_candidate=max(rows, key=lambda r: r['candidate_correct']),
                              best_unique=max(decisions, key=lambda r: (r['correct'], -r['errors'])),
                              best_budget=min(decisions, key=lambda r: (r['max_lane_budget'], -r['correct'], r['errors']))))
        status = Path(summaries[-1]['best_candidate']['candidate_path'].replace('.truth_candidates.tsv', '.status.tsv'))
        with status.open() as f:
            next(f)
            summaries[-1]['best_candidate_source_status_counts'] = dict(Counter(row.rstrip().split('\t')[-1] for row in itertools.islice(f, START, N)))
    summary = dict(plan=plan, complete=completed == expected, completed_new_configurations=len(completed),
                   baseline_candidate_union=nbits(baseline), new_candidate_union=nbits(new), combined_candidate_union=nbits(combined),
                   combined_recall=nbits(combined)/B, newly_recovered=nbits(new & ~baseline), still_missing=B-nbits(combined),
                   zero_candidate_union=nbits(zero), only_inner_rs=nbits(combined & ~zero),
                   best_single_candidate=best_candidate, best_single_unique=best_unique, best_outer_budget=best_budget,
                   baseline_best_candidate=old_group['best_single_candidate'], baseline_best_unique=old_group['best_single_unique'],
                   baseline_best_outer_budget=old_group['best_outer_budget'], by_k=summaries,
                   interpretation='Truth-labelled index-wise unions are evaluation upper coverage. They do not constitute one reference-blind decoder output. Read assignment used the full design.',
                   metrics='Each correct candidate equals the entire designed 60-nt strand after inner RS. One index counts once across all configurations.')
    if args.partial:
        (OUT / 'partial_summary.json').write_text(json.dumps(summary, indent=2))
        print('PARTIAL', len(completed), 'best', best_candidate, 'union', nbits(combined), flush=True)
        return
    for name, rows in [('parameters.tsv', configs), ('decision_metrics.tsv', scores)]:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (OUT / name).open('w') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    matrix = np.stack(matrices, axis=1)
    before = np.any(matrix[:, [j for j, r in enumerate(configs) if r['experiment'] == 'baseline']], axis=1)
    after = np.any(matrix, axis=1)
    assert before.sum() == nbits(baseline) and after.sum() == nbits(combined)
    with (OUT / 'per_index.tsv').open('w') as f:
        f.write('index\tbaseline179_success\tcombined_success\tnewly_recovered\tsuccessful_configurations\n')
        for j in range(B): f.write(f'{START+j}\t{int(before[j])}\t{int(after[j])}\t{int(after[j] and not before[j])}\t{int(matrix[j].sum())}\n')
    with gzip.open(OUT / 'candidate_success_matrix.tsv.gz', 'wt') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['index'] + [f"{r['experiment']}:k{r['k']}.c{r['cutoff']}.r{r['ratio']:.2f}" for r in configs])
        writer.writerows([START+j, *matrix[j]] for j in range(B))
    # Export the actual lowest-error-budget output, retaining errors and filtering only by numeric block.
    selected = OUT / 'best_budget.strands'; fasta = OUT / 'best_budget.fasta'
    with Path(best_budget['path']).open() as src, selected.open('w') as dst, fasta.open('w') as fa:
        for line in src:
            field, dna = line.rstrip().split('\t'); index = int(field)
            if START <= index < N: dst.write(line); fa.write(f'>{index}\n{dna}\n')
    checked, _ = score_outputs(selected, truth, [12])
    for field in ['correct', 'errors', 'erasures', 'lane_errors', 'lane_budgets']:
        assert checked[12][field] == best_budget[field]
    control_dir = OUT / 'reference_blind_control'; control_dir.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update(ANTK_LAYOUT='internal', ANTK_BLOOM='1', ANTK_READ_LIMIT='200',
               ANTK_GRID=f"{best_budget['ratio']:.2f}:{best_budget['cutoff']}")
    command = [str(DATA / 'bin/graph_block_rs'), plan['input'], str(N), str(best_budget['k']), str(control_dir / 'file3'), '-']
    with (control_dir / 'run.log').open('w') as log:
        subprocess.run(command, env=env, stdout=log, stderr=log, check=True)
    label = f"file3.k{best_budget['k']}.c{best_budget['cutoff']}.r{best_budget['ratio']:.2f}"
    original_stem = Path(best_budget['path']).parent / label
    control_stem = control_dir / label
    policy_names = json.loads(Path(str(original_stem) + '.inner_stats.json').read_text())['emitted']
    checks = {}
    for suffix in [f'.{p}.strands' for p in policy_names] + ['.status.tsv', '.decisions.tsv', '.inner_stats.json']:
        original = Path(str(original_stem) + suffix).read_bytes()
        assert original == Path(str(control_stem) + suffix).read_bytes(), suffix
        checks[suffix] = hashlib.sha256(original).hexdigest()
    summary['reference_blind_control'] = dict(command=command, identical_sha256=checks, passed=True,
        scope='Graph traversal and inner-RS decisions receive no design sequences. Earlier read assignment remains reference-assisted.')
    summary['selected_output'] = str(selected); summary['selected_fasta'] = str(fasta)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2))
    plot(summary, configs, scores)
    lines = ['# 第 12 块小 k-mer 扩展测试', '',
             f"同一批参考辅助分块的 {plan['input_reads']:,} 条原始 reads；长度上限 200、路径上限 1050、全局 index 域 196596、内码 RS 与历史实验一致。`-c c` 的有效计数条件为 count > c。", '',
             f"新增 {len(completed)} 组参数，与历史 179 组合计 {len(configs)} 组。", '',
             '| k | 合计参数组 | 历史候选并集 | 扩展后同 k 并集 | 最佳单参数正确候选 | 相对历史全部参数新增 |',
             '|---|---:|---:|---:|---:|---:|']
    for g in summaries:
        lines.append(f"| {g['k']} | {g['parameters']} | {g['baseline_union']} | {g['combined_union']} | {g['best_candidate']['candidate_correct']} | {g['newly_recovered_vs_all_baseline']} |")
    lines += ['', f"全部参数并集：{nbits(baseline)} → **{nbits(combined)} / {B} ({nbits(combined)/B:.4%})**，新增 {nbits(new & ~baseline)}，仍缺 {B-nbits(combined)}。其中仅通过内码 RS 才进入全局并集的片段为 {nbits(combined & ~zero)}。", '',
              f"若将并集中的完整正确片段全部选出、其余均作为擦除，擦除数仍为 {B-nbits(combined)}；纯擦除恢复要求至少 10977 条正确片段，当前距离该条件为 {max(0, 10977-nbits(combined))} 条。这一计算不使用错误片段中可能仍正确的部分 lane。", '',
              f"最佳单参数候选：k={best_candidate['k']}, c={best_candidate['cutoff']}, r={best_candidate['ratio']}，正确候选 {best_candidate['candidate_correct']}。", '',
              f"最佳实际唯一输出：k={best_unique['k']}, c={best_unique['cutoff']}, r={best_unique['ratio']}, {best_unique['policy']}；正确 {best_unique['correct']}，错误 {best_unique['errors']}，擦除 {best_unique['erasures']}。", '',
              f"最低外码预算实际输出：k={best_budget['k']}, c={best_budget['cutoff']}, r={best_budget['ratio']}, {best_budget['policy']}；正确 {best_budget['correct']}，错误 {best_budget['errors']}，擦除 {best_budget['erasures']}。", '',
              f"6 路错误数：{best_budget['lane_errors']}；2E+S：{best_budget['lane_budgets']}；最差 {best_budget['max_lane_budget']}，要求 ≤5406。", '',
              '将最优外码预算参数的设计真值参数替换为 `-` 再运行，全部 11 个实际输出策略、status、decisions、inner_stats 文件均逐字节一致，验证图搜索及内码 RS 决策不依赖评分真值。该对照不改变分块过程使用设计序列的事实。', '',
              '候选并集通过设计真值判断每个 index 在任一参数下是否出现完整正确的 60 nt RS 码字，仅作可恢复覆盖上界；不能当作一个无需真值的实际解码输出。实际输出由图与内码 RS 独立选择后再评分，导出的 best_budget 文件保留了所有错误。分块本身使用了完整设计序列，属于参考辅助实验。', '',
              '逐 index 的参数成功情况：candidate_success_matrix.tsv.gz（1 表示完整正确候选存在，0 表示不存在），汇总见 per_index.tsv；逐参数统计见 parameters.tsv；所有实际决策策略的正确/错误与外码预算见 decision_metrics.tsv。', '']
    (OUT / 'REPORT.zh-CN.md').write_text('\n'.join(lines))
    names = ['summary.json', 'parameters.tsv', 'decision_metrics.tsv', 'per_index.tsv', 'candidate_success_matrix.tsv.gz', 'best_budget.strands', 'best_budget.fasta', 'REPORT.zh-CN.md', 'comparison.png', 'comparison.svg']
    (OUT / 'validation.json').write_text(json.dumps(dict(candidate_matrix_shape=list(matrix.shape), baseline_union_verified=10033,
        actual_output_rescored=checked[12], source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), Path(__file__).with_name('run_block12_small_k.py')]},
        artifact_sha256={name: hashlib.sha256((OUT / name).read_bytes()).hexdigest() for name in names}), indent=2))
    print(json.dumps({k: summary[k] for k in ['completed_new_configurations', 'combined_candidate_union', 'newly_recovered', 'best_single_candidate', 'best_outer_budget']}, indent=2), flush=True)


def plot(summary, configs, scores):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    groups = [r for r in summary['by_k'] if r['k'] <= 21]
    x = [r['k'] for r in groups]
    axes[0].plot(x, [r['combined_union'] for r in groups], 'o-', label='Union across parameters at this k')
    axes[0].plot(x, [r['best_candidate']['candidate_correct'] for r in groups], 's-', label='Best single-parameter candidate coverage')
    axes[0].plot(x, [r['best_unique']['correct'] for r in groups], '^--', label='Best actual unique-output correct count')
    axes[0].set(xlabel='k-mer size', ylabel='Correct strands / 16,383', title='Block 12: smaller k has an optimum', xticks=x)
    rows = sorted([r for r in configs if r['k'] == 16 and r['cutoff'] == 1 and r['ratio'] <= 5], key=lambda r: r['ratio'])
    axes[1].plot([r['ratio'] for r in rows], [r['candidate_correct'] for r in rows], 'o-', label='Correct candidates after inner RS')
    strict = sorted([r for r in scores if r['k'] == 16 and r['cutoff'] == 1 and r['ratio'] <= 5 and r['policy'] == 'strict'], key=lambda r: r['ratio'])
    axes[1].plot([r['ratio'] for r in strict], [r['correct'] for r in strict], 's--', label='Correct actual strict output')
    axes[1].set(xlabel='Adjacent k-mer coverage ratio (-r)', ylabel='Correct strands / 16,383', title='k=16, cutoff=1: ratio refinement')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.13), fontsize=8, frameon=False)
    fig.suptitle('Same 1,681,943 reference-assigned reads; identical search limits', fontsize=11)
    fig.savefig(OUT / 'comparison.png', dpi=180)
    fig.savefig(OUT / 'comparison.svg')
    plt.close(fig)


if __name__ == '__main__': main()
