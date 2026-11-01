#!/usr/bin/env python3
"""Run the registered k / cutoff / ratio grid and collect reference-only metrics."""
import argparse
import concurrent.futures as cf
import csv
import json
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / 'input/antkowiak2020'


def truth_sequences(path):
    return [s.strip() for s in path.read_text().splitlines()
            if s.strip() and not s.startswith('>')]


def summarize(label):
    truth_path=DATA/f'sources/File{label[-1]}_design.fasta'
    if label!='file1':truth_path=truth_path.with_name(f'File{label[-1]}_synthesized.fasta')
    truth = truth_sequences(truth_path)
    results = []
    support_cache = {}
    for p in sorted((DATA/'results').glob(label + '.k*.strands')):
        status = list(csv.DictReader(p.with_suffix('.status.tsv').open(), delimiter='\t'))
        stem = p.name.split('.')
        k = int(stem[1][1:])
        support_path=DATA/f'results/{label}.k{k}.truth_support.tsv'
        if len(status)!=len(truth) or not support_path.exists():
            print('Skipping unfinished search/diagnostic:',p.name,file=sys.stderr)
            continue
        if k not in support_cache:
            with support_path.open() as handle:
                support_cache[k]=list(csv.DictReader(handle,delimiter='\t'))
        support=support_cache[k]
        if len(support)!=len(truth):
            print('Skipping unfinished support diagnostic:',p.name,file=sys.stderr)
            continue
        output = dict((int(line.split('\t')[0]), line.split('\t')[1])
                      for line in p.read_text().splitlines())
        correct = sum(0 <= i < len(truth) and s == truth[i] for i, s in output.items())
        # Ratio contains a dot, so read it from the remaining name.
        ratio = float(p.name.split('.r', 1)[1].removesuffix('.strands'))
        cutoff = int(stem[2][1:])
        block_errors = []
        for start in range(0, len(truth), 16383):
            count = min(16383, len(truth)-start)
            emitted = sum(start <= i < start+count for i in output)
            errors = sum(start <= i < start+count and s != truth[i] for i, s in output.items())
            block_errors.append(dict(block=start//16383, errors=errors,
                                     erasures=count-emitted, budget=2*errors+count-emitted,
                                     within_rs_bound=2*errors+count-emitted <= 5406))
        row = dict(dataset=label, k=k, cutoff=cutoff, ratio=ratio,
                   truth=len(truth), emitted=len(output), correct=correct,
                   errors=len(output)-correct, erasures=len(truth)-len(output),
                   exact_recall=correct/len(truth), precision=correct/len(output) if output else None,
                   overflow=sum(s['status']=='overflow' for s in status),
                   ambiguous=sum(s['status']=='ambiguous' for s in status),
                   graph_supported=sum(int(s['min_kmer_count'])>cutoff for s in support),
                   per_block=block_errors)
        results.append(row)
    (DATA/f'results/{label}_graph_metrics.json').write_text(json.dumps(results, indent=2))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dataset', choices=['file1', 'file2', 'file3'])
    ap.add_argument('--k', nargs='+', type=int, default=[12, 13, 15, 17, 21])
    ap.add_argument('--parallel', type=int, default=3)
    ap.add_argument('--summarize-only', action='store_true')
    ap.add_argument('--cutoffs')
    ap.add_argument('--ratio', default='2.55', help='One ratio or comma-separated ratios; reuses each k graph')
    args = ap.parse_args()
    if args.summarize_only:
        summarize(args.dataset)
        return
    label = args.dataset
    ready = DATA/f'results/{label}_reads.json'
    while not ready.exists() or not ready.stat().st_size:
        time.sleep(5)
    truth = DATA/f'sources/File{label[-1]}_design.fasta'
    if label!='file1':truth=truth.with_name(f'File{label[-1]}_synthesized.fasta')
    n = len(truth_sequences(truth))

    def run(k):
        prefix = DATA/f'results/{label}'
        command = [str(DATA/'bin/graph_adapter'), str(DATA/f'processed/{label}.fa.gz'),
                   str(n), str(k), str(prefix), str(truth)]
        if args.cutoffs:
            command += [args.ratio, args.cutoffs]
        env = os.environ.copy()
        env['ANTK_LAYOUT'] = 'internal'
        env['ANTK_BLOOM'] = '0' if label=='file1' else '1'
        start = time.monotonic()
        logname = DATA/f'results/{label}_k{k}'
        if args.cutoffs:
            logname = pathlib.Path(str(logname) + '_cutoffs_' + args.cutoffs.replace(',', '-'))
            logname = pathlib.Path(str(logname) + '_ratios_' + args.ratio.replace(',', '-'))
        with pathlib.Path(str(logname)+'.log').open('w') as actual_log:
            p = subprocess.run(['/usr/bin/time', '-l']+command, stdout=actual_log, stderr=actual_log, env=env)
        result = dict(command=command, k=k, exit_code=p.returncode, seconds=time.monotonic()-start)
        pathlib.Path(str(logname)+'.run.json').write_text(json.dumps(result, indent=2))
        print('FINISHED', label, json.dumps(result), flush=True)
        if p.returncode:
            raise RuntimeError(result)

    with cf.ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(run, args.k))
    print(json.dumps(summarize(label), indent=2), flush=True)


if __name__ == '__main__':
    main()
