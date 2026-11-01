#!/usr/bin/env python3
"""Summarize the completed File 3 ratio sweep, with the original runs retained."""
import csv
import json
import sys

from report import DATA, OUT, load, truth
from run_graphs import summarize


def rank(row):
    return row['correct'], -row['errors'], -row['overflow']


def filename(row):
    return OUT / (f"file3.k{row['k']}.c{row['cutoff']}.r{row['ratio']:.2f}.strands")


def exact_ids(row, target):
    return {int(i) for i, dna in (line.split('\t') for line in filename(row).read_text().splitlines())
            if target[int(i)] == dna}


def main():
    rows = summarize('file3')
    target = truth('file3')
    # A status file must contain exactly one decision per design index. Avoid
    # accidentally reporting partial outputs from a still-running search.
    for row in rows:
        with filename(row).with_suffix('.status.tsv').open() as handle:
            assert sum(1 for _ in handle) == len(target) + 1, filename(row)
    old_rows = load(OUT / 'file3_before_ratio_sweep.json')
    grids=[([15,17,19,21],[1,3,9,19],[1.2,1.5,1.8,2.1,3.5]),
           ([13],[9,49,99],[1.2,1.5,1.8]),
           ([18,20],[1,3,9],[1.8,2.1,2.55,3,3.5,5,10]),
           ([19],[1,3,9],[3,5,10]),
           ([18],[2],[1.8,2.1,2.55,3.5])]
    expected={(k,c,r) for ks,cs,rs in grids for k in ks for c in cs for r in rs}
    present={(r['k'],r['cutoff'],r['ratio']) for r in rows}
    assert expected<=present, ('Follow-up configurations not complete',sorted(expected-present))
    old = max(old_rows, key=rank)
    best = max(rows, key=rank)
    block_best=max(rows,key=lambda r:(sum(b['within_rs_bound'] for b in r['per_block']),)+rank(r))
    by_k = [max((r for r in rows if r['k'] == k), key=rank)
            for k in sorted({r['k'] for r in rows})]
    by_ratio = [max((r for r in rows if r['ratio'] == ratio), key=rank)
                for ratio in sorted({r['ratio'] for r in rows})]
    before, after = exact_ids(old, target), exact_ids(best, target)
    ceilings=[]
    for k in sorted({r['k'] for r in rows}):
        with (OUT/f'file3.k{k}.truth_support.tsv').open() as handle:
            support=list(csv.DictReader(handle,delimiter='\t'))
        block_counts=[sum(int(s['min_kmer_count'])>1 for s in support[start:start+16383])
                      for start in range(0,len(target),16383)]
        ceilings.append(dict(k=k,cutoff=1,per_block_true_path_count=block_counts,
                             weakest_block_count=min(block_counts),
                             last_block_count=block_counts[-1],
                             necessary_all_blocks_supported=all(n>=10977 for n in block_counts)))
    record = dict(configurations=len(rows), new_configurations=len(rows)-len(old_rows),
                  followup_grids=[dict(k=ks,cutoffs=cs,ratios=rs) for ks,cs,rs in grids],
                  original_best=old, best=best, best_block_feasibility=block_best, best_per_k=by_k,
                  best_per_ratio=by_ratio, gained_exact=len(after-before),
                  lost_exact=len(before-after),
                  recall_gain_percentage_points=100*(best['exact_recall']-old['exact_recall']),
                  graph_support_ceilings=ceilings,
                  all_blocks_in_bound=sum(all(b['within_rs_bound'] for b in r['per_block']) for r in rows))
    (OUT / 'file3_ratio_summary.json').write_text(json.dumps(record, indent=2))
    fields=['k','cutoff','ratio','correct','errors','erasures','exact_recall','overflow','ambiguous']
    with (OUT / 'file3_ratio_comparison.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields+['blocks_in_bound'],extrasaction='ignore')
        writer.writeheader()
        writer.writerows(dict(r,blocks_in_bound=sum(b['within_rs_bound'] for b in r['per_block'])) for r in rows)
    sys.path.insert(0,str(DATA/'plotlibs'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax=plt.subplots(figsize=(8.2,4.7),layout='constrained')
    for k in sorted({r['k'] for r in rows}):
        points=[max((r for r in rows if r['k']==k and r['ratio']==ratio),key=rank)
                for ratio in sorted({r['ratio'] for r in rows if r['k']==k})]
        ax.plot([p['ratio'] for p in points],[100*p['exact_recall'] for p in points],'-o',label=f'k={k}')
    ax.set(xlabel='Coverage ratio (-r)',ylabel='Exact strand recall (%)',
           title='File 3: best tested cutoff at each k and ratio',ylim=(-2,100))
    ax.set_xscale('log')
    ticks=sorted({r['ratio'] for r in rows})
    ax.set_xticks(ticks,[f'{r:g}' for r in ticks],rotation=45)
    ax.minorticks_off()
    ax.grid(alpha=.2)
    ax.legend(frameon=False,ncol=2)
    fig.savefig(OUT/'file3_ratio_comparison.png',dpi=180)
    fig.savefig(OUT/'file3_ratio_comparison.svg')
    plt.close(fig)
    text=['**File 3：小 k 与覆盖比参数补充实验**','',
          f"原实验 {len(old_rows)} 组，本次新增 {record['new_configurations']} 组；合计 {len(rows)} 组。使用全部 92,337,721 条保留 reads，196,596 条参考链。程序仍为实验适配版，路径上限 1050，严格内层 RS 校验；参考仅在搜索后用于评分。",'',
          '`-r` 控制相邻 k-mer 覆盖数允许的比值；越小通常越严格。沿用 DBGPS2 的整数截断和饱和规则，故并非精确的实数对称比值条件。覆盖截断 c 要求计数严格大于 c。每次运行先对一个 k 计数，再复用该图扫描多个 ratio 和 c；没有改变图搜索及内码判定。','',
          '| 范围 | k | c | -r | 正确链 | 恢复率 | 错误链 | 溢出 | 满足外码界的块 |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for label,r in [('原最佳',old),('当前最佳',best)]:
        text.append(f"| {label} | {r['k']} | {r['cutoff']} | {r['ratio']:.2f} | {r['correct']:,} | {100*r['exact_recall']:.2f}% | {r['errors']} | {r['overflow']:,} | {sum(b['within_rs_bound'] for b in r['per_block'])}/12 |")
    text += ['',f"按外码块达标数量选优，当前最多为 {sum(b['within_rs_bound'] for b in block_best['per_block'])}/12 块（k={block_best['k']}、c={block_best['cutoff']}、ratio={block_best['ratio']:.2f}）。"]
    text += ['',f"相对原最佳，正确恢复率变化 {record['recall_gain_percentage_points']:+.2f} 个百分点；新得到 {record['gained_exact']:,} 条正确链，原来正确但当前最佳未正确输出 {record['lost_exact']:,} 条。这里比较单一参数配置，没有合并不同配置的输出。",'',
             '**每个 k 的已测最佳结果**','',
             '| k | c | -r | 正确链 | 恢复率 | 错误链 | 溢出 |',
             '|---:|---:|---:|---:|---:|---:|---:|']
    for r in by_k:
        text.append(f"| {r['k']} | {r['cutoff']} | {r['ratio']:.2f} | {r['correct']:,} | {100*r['exact_recall']:.2f}% | {r['errors']} | {r['overflow']:,} |")
    text += ['', '**固定 k=19、c=1，仅改变 ratio**','',
             '| -r | 正确链 | 恢复率 | 错误链 | 溢出 | 满足外码界的块 |',
             '|---:|---:|---:|---:|---:|---:|']
    for r in sorted((r for r in rows if r['k']==19 and r['cutoff']==1),key=lambda r:r['ratio']):
        text.append(f"| {r['ratio']:.2f} | {r['correct']:,} | {100*r['exact_recall']:.2f}% | {r['errors']} | {r['overflow']:,} | {sum(b['within_rs_bound'] for b in r['per_block'])}/12 |")
    text += ['',f"![File 3 参数扫描]({OUT/'file3_ratio_comparison.png'})",'',
             '曲线在每个 k/r 组合中选择已测最佳 c，各组合的 c 网格可能不同；未测组合不计入结果。所有“最佳”均依据本数据真值事后选择，属于探索结果，没有独立验证集。','',
             '小 k 会丢失序列辨识能力：k=13 的计数日志记录了 786382 个计数大于 1 的索引锚点。在 196596 个已知索引、12 nt 索引长度下，可能的 13-mer 锚点总数仅为 196596×4=786384，几乎全部被观察到。这种锚点上下文的饱和与完整正确链是否可唯一恢复是两回事。','',
             '**当前最佳的逐块外码界**','',
             '| 块（从 1 起） | 错误 E | 缺失 S | 2E+S | ≤5406 |',
             '|---:|---:|---:|---:|---|']
    for b in best['per_block']:
        text.append(f"| {b['block']+1} | {b['errors']} | {b['erasures']} | {b['budget']} | {'是' if b['within_rs_bound'] else '否'} |")
    text += ['', '**仅调 ratio 无法补回图中不存在的真值路径**','',
             '以下在最宽松的已测覆盖截断 c=1 下，统计全部真实 k-mer 计数均大于 1 的参考链数，不施加 ratio 或路径上限。它是当前严格全路径搜索能正确输出的上限，不是实际恢复数。每个块至少需要 10977 条正确链（有错误时要求更多），才能满足这里采用的标准 RS 错误/缺失界。','',
             '| k | 第 12 块的完整真值路径数 | 所有块中最小值 | 每块均达到 10977 |',
             '|---:|---:|---:|---|']
    for c in ceilings:
        text.append(f"| {c['k']} | {c['last_block_count']:,} | {c['weakest_block_count']:,} | {'是，仅满足必要条件' if c['necessary_all_blocks_supported'] else '否'} |")
    conclusion='本扫描未取得完整文件恢复' if not record['all_blocks_in_bound'] else '满足纠错界后仍需实际解码并验证文件'
    text += ['',f"全部已测配置中，有 {record['all_blocks_in_bound']} 组的 12 个块均满足标准 RS 错误/缺失纠错界。{conclusion}；参考序列解出的 ZIP 不属于测序恢复结果。",'',
             f"完整参数与指标：[CSV]({OUT/'file3_ratio_comparison.csv'})；汇总：[JSON]({OUT/'file3_ratio_summary.json'})。逐索引输出、状态、计数及运行日志保存在同一目录。未修改生产 DBGPS2 源码，未提交或推送。"]
    (OUT/'FILE3_RATIO_REPORT.zh-CN.md').write_text('\n'.join(text)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ['best_per_k','best_per_ratio']},indent=2))


if __name__=='__main__':
    main()
