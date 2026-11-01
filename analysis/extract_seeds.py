#!/usr/bin/env python3
"""Extract droplet seed observations from the (already adapter-trimmed)
DNA Fountain FASTQ files deposited at ENA.

R1 starts exactly at the droplet block: seed(16nt) + payload + RS (151 nt).
R2 is the mate; its reverse complement starts one base into the block, so
only R1 yields direct seed reads. Writes `seed\tcount`, descending count.
"""
import argparse
import gzip
import sys
from collections import Counter

from fountain_common import dna_to_bytes

LFSR_POLY = 0b100000000000000000000000011000101
LFSR_START = 0b001010101


def lfsr_states(n):
    """First n states of the encoder's 32-bit Galois LFSR, in order."""
    nbits = LFSR_POLY.bit_length() - 1
    state = LFSR_START
    out = []
    for _ in range(n):
        state <<= 1
        if (state >> nbits) & 1:
            state ^= LFSR_POLY
        out.append(state)
    return out


def build_seed_set(n=400000):
    return set(lfsr_states(n))


def scan_r1(path, counter, valid):
    op = gzip.open if path.endswith('.gz') else open
    n = 0
    with op(path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            block = line.strip()[:16]
            if len(block) == 16 and set(block) <= set('ACGT'):
                seed = int.from_bytes(dna_to_bytes(block), 'big')
                if seed in valid:
                    counter[seed] += 1
                    n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--r1', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--min_count', type=int, default=1,
                    help='only emit seeds observed at least this often')
    ap.add_argument('--lfsr_states', type=int, default=600000,
                    help='number of leading LFSR states accepted as valid seeds')
    args = ap.parse_args()
    valid = build_seed_set(args.lfsr_states)
    counter = Counter()
    n = scan_r1(args.r1, counter, valid)
    with open(args.out, 'w') as f:
        for seed, cnt in counter.most_common():
            if cnt >= args.min_count:
                f.write(f"{seed}\t{cnt}\n")
    print(f"R1 seed observations: {n}; distinct seeds: {len(counter)}; "
          f"kept (>={args.min_count}): {sum(1 for c in counter.values() if c >= args.min_count)}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
