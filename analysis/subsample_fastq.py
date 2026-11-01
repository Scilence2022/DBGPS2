#!/usr/bin/env python3
"""Subsample paired FASTQ.gz by read pairs (random, reproducible seed),
writing R1/R2 subsets as plain gzip FASTQ for DBGPS2 input."""
import argparse
import gzip
import random


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--r1', required=True)
    ap.add_argument('--r2', required=True)
    ap.add_argument('--pairs', type=int, required=True, help='number of read pairs to keep')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--out1', required=True)
    ap.add_argument('--out2', required=True)
    args = ap.parse_args()

    # Count pairs first (cheap line counting)
    def count_records(path):
        n = 0
        with gzip.open(path, 'rb') as f:
            for i, _ in enumerate(f):
                pass
        assert (i + 1) % 4 == 0
        return (i + 1) // 4

    total = min(count_records(args.r1), count_records(args.r2))
    rng = random.Random(args.seed)
    keep = set(rng.sample(range(total), min(args.pairs, total)))
    kept = 0
    with gzip.open(args.r1, 'rt') as f1, gzip.open(args.r2, 'rt') as f2, \
            gzip.open(args.out1, 'wt') as o1, gzip.open(args.out2, 'wt') as o2:
        for idx in range(total):
            r1 = [f1.readline() for _ in range(4)]
            r2 = [f2.readline() for _ in range(4)]
            if idx in keep:
                o1.writelines(r1)
                o2.writelines(r2)
                kept += 1
    print(f"kept {kept} pairs of {total}")


if __name__ == '__main__':
    main()
