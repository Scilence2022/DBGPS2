#!/usr/bin/env python3
"""Paper-style prefix clustering, MUSCLE 3.8.31, and weighted voting.

Same reads as graph experiments; prefix=16, cluster minimum=5, up to 15 reads,
five uniform samples for large clusters, gap vote weight=0.4. This follows the
article's parameters (the released notebook instead uses 14/4/0.5).
"""
import argparse
import collections
import concurrent.futures as cf
import gzip
import json
import pathlib
import random
import subprocess
import time
COMP = str.maketrans('ACGT', 'TGCA')
def canonical(s):
    return min(s, s.translate(COMP)[::-1])


def process(job):
    if len(job) == 3:
        key, reads, muscle = job
        gap_weight = .4
    else:
        key, reads, muscle, gap_weight = job
    rng = random.Random(20261017 ^ int.from_bytes(key.encode(), 'big'))
    samples = [reads] if len(reads) <= 15 else [rng.sample(reads, 15) for _ in range(5)]
    candidates = []
    for sample in samples:
        fasta = ''.join(f'>s{i}\n{s}\n' for i, s in enumerate(sample))
        p = subprocess.run([muscle, '-quiet'], input=fasta,
                           capture_output=True, text=True, timeout=60, check=True, close_fds=False)
        aligned = []
        for line in p.stdout.splitlines():
            if line.startswith('>'):
                aligned.append('')
            else:
                aligned[-1] += line.strip().upper()
        assert len(aligned) == len(sample) and len(set(map(len, aligned))) == 1
        sequence = []
        for col in zip(*aligned):
            counts = dict.fromkeys('ACGT-N', 0)
            for c in col:
                counts[c] += 1
            counts['-'] *= gap_weight
            c = max(counts, key=counts.get)
            if c != '-':
                sequence.append(c)
        seq = ''.join(sequence)
        candidates.append(seq[:60] if len(seq) >= 60 else seq)
    return candidates


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('reads')
    ap.add_argument('truth')
    ap.add_argument('muscle')
    ap.add_argument('prefix')
    ap.add_argument('--workers', type=int, default=12)
    ap.add_argument('--prefix-length', type=int, default=16)
    ap.add_argument('--min-cluster', type=int, default=5)
    ap.add_argument('--gap-weight', type=float, default=.4)
    ap.add_argument('--raw-fastq', action='store_true')
    args = ap.parse_args()
    truth = {canonical(s.strip()) for s in pathlib.Path(args.truth).read_text().splitlines()
             if s.strip() and not s.startswith('>')}
    start = time.monotonic()
    groups = collections.defaultdict(list)
    raw_exact = set()
    total = 0
    with gzip.open(args.reads, 'rt') as f:
        for line in f:
            if args.raw_fastq:
                seq = next(f).strip()
                next(f)
                next(f)
                if not 55 <= len(seq) <= 70:
                    continue
            else:
                if line.startswith('>'):
                    continue
                seq = line.strip()
            if canonical(seq[:60]) in truth:
                raw_exact.add(canonical(seq[:60]))
            groups[seq[:args.prefix_length]].append(seq)
            total += 1
    jobs = [(key, value) for key, value in groups.items() if len(value) >= args.min_cluster]
    jobs.sort()
    del groups
    cluster_seconds = time.monotonic() - start
    print('CLUSTERED', total, len(jobs), 'raw exact', len(raw_exact),
          'seconds', round(cluster_seconds, 3), flush=True)


    recovered = set()
    generated = written = exact = 0
    msa_start = time.monotonic()
    with open(args.prefix + '.candidates.txt', 'w') as out, cf.ProcessPoolExecutor(args.workers) as pool:
        for done, candidates in enumerate(pool.map(process, ((k, v, args.muscle, args.gap_weight) for k, v in jobs), chunksize=64), 1):
            for seq in candidates:
                generated += 1
                if len(seq) == 60 and set(seq) <= set('ACGT'):
                    written += 1
                    out.write(seq + '\n')
                    if canonical(seq) in truth:
                        exact += 1
                        recovered.add(canonical(seq))
            if done % 5000 == 0:
                print('MSA', done, '/', len(jobs), 'recovered', len(recovered),
                      'seconds', round(time.monotonic() - msa_start, 1), flush=True)
    result = dict(method=f'prefix{args.prefix_length}_muscle3_vote{args.gap_weight}', reads=total, truth=len(truth),
                  prefix_length=args.prefix_length, min_cluster=args.min_cluster,
                  gap_weight=args.gap_weight, raw_fastq=args.raw_fastq,
                  clusters=len(jobs), raw_exact_recovered=len(raw_exact),
                  candidates_generated=generated, candidates_60nt=written,
                  exact_candidates=exact, recovered=len(recovered),
                  recall=len(recovered)/len(truth), seed=20261017,
                  cluster_seconds=cluster_seconds,
                  msa_seconds=time.monotonic()-msa_start,
                  total_seconds=time.monotonic()-start)
    pathlib.Path(args.prefix + '.json').write_text(json.dumps(result, indent=2))
    print('RESULT', json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
