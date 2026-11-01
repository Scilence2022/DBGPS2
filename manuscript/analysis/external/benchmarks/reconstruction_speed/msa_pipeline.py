"""Reference-free prefix clustering, MUSCLE 3.8.31, weighted majority vote.

Uses the previously reported photolithographic baseline parameters. No design
sequence is loaded. Timed output includes inner-code selection; truth scoring
is performed afterwards by run_benchmark.py.
"""
import argparse
import collections
import concurrent.futures
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import time


def reconstruct(job):
    key, reads, muscle = job
    rng = random.Random(20261017 ^ int.from_bytes(key.encode(), 'big'))
    subsets = [reads] if len(reads) <= 15 else [rng.sample(reads, 15) for _ in range(5)]
    candidates = []
    alignment_cpu = voting_cpu = 0.0
    for subset in subsets:
        fasta = ''.join(f'>s{i}\n{s}\n' for i, s in enumerate(subset))
        start = time.perf_counter()
        p = subprocess.run([muscle, '-quiet'], input=fasta, capture_output=True,
                           text=True, timeout=60, check=True, close_fds=False)
        alignment_cpu += time.perf_counter() - start
        start = time.perf_counter()
        aligned = []
        for line in p.stdout.splitlines():
            if line.startswith('>'):
                aligned.append('')
            else:
                aligned[-1] += line.strip().upper()
        assert len(aligned) == len(subset)
        assert len(set(map(len, aligned))) == 1
        sequence = []
        for column in zip(*aligned):
            counts = dict.fromkeys('ACGT-N', 0)
            for base in column:
                counts[base] += 1
            counts['-'] *= .4
            base = max(counts, key=counts.get)
            if base != '-':
                sequence.append(base)
        sequence = ''.join(sequence)
        if len(sequence) >= 60:
            sequence = sequence[:60]
        if len(sequence) == 60 and set(sequence) <= set('ACGT'):
            candidates.append(sequence)
        voting_cpu += time.perf_counter() - start
    return candidates, len(subsets), alignment_cpu, voting_cpu


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reads')
    parser.add_argument('design_count', type=int)
    parser.add_argument('output_prefix')
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    prefix = Path(args.output_prefix)
    start = time.perf_counter()
    groups = collections.defaultdict(list)
    reads = 0
    with gzip.open(args.reads, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                continue
            sequence = line.strip()
            groups[sequence[:16]].append(sequence)
            reads += 1
    jobs = sorted((k, v) for k, v in groups.items() if len(v) >= 5)
    group_count = len(groups)
    del groups
    clustered = time.perf_counter()
    print(f'CLUSTERED reads={reads} eligible={len(jobs)} seconds={clustered-start:.3f}', flush=True)
    candidates = alignments = 0
    alignment_worker_seconds = voting_worker_seconds = 0.0
    with open(str(prefix)+'.candidates.txt', 'w') as output:
        with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
            iterator = pool.map(reconstruct,
                                ((k, v, str(root/'bin/muscle3')) for k, v in jobs),
                                chunksize=64)
            for done, (sequences, count, align_time, vote_time) in enumerate(iterator, 1):
                for sequence in sequences:
                    output.write(sequence+'\n')
                candidates += len(sequences)
                alignments += count
                alignment_worker_seconds += align_time
                voting_worker_seconds += vote_time
                if done % 25000 == 0:
                    print(f'MSA done={done}/{len(jobs)} seconds={time.perf_counter()-clustered:.1f}', flush=True)
    voted = time.perf_counter()
    selection = subprocess.run([str(root/'bin/inner_select'),
                                str(prefix)+'.candidates.txt', str(args.design_count),
                                str(prefix)], check=True, capture_output=True, text=True)
    finished = time.perf_counter()
    record = dict(reads=reads, design_count=args.design_count, workers=args.workers,
                  prefix_length=16, min_cluster=5, max_alignment_reads=15,
                  large_cluster_subsets=5, gap_weight=.4, seed=20261017,
                  prefix_groups=group_count, eligible_clusters=len(jobs),
                  alignment_calls=alignments, candidates_60nt=candidates,
                  cluster_seconds=clustered-start, msa_vote_seconds=voted-clustered,
                  validation_seconds=finished-voted, pipeline_seconds=finished-start,
                  alignment_worker_elapsed_sum=alignment_worker_seconds,
                  voting_worker_elapsed_sum=voting_worker_seconds,
                  selection=json.loads(selection.stdout), reference_access=False)
    Path(str(prefix)+'.stages.json').write_text(json.dumps(record, indent=2))
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()
