#!/usr/bin/env python3
"""Position-wise consensus of the droplet block per seed from R1 reads
(which cover block[0:151]; the last RS base is left to RS validation).
Key: full 16-nt seed. Output: seed, consensus, depth, min position depth."""
import argparse
import gzip
import sys
from collections import Counter, defaultdict

from fountain_common import dna_to_bytes

BLOCK_FROM_R1 = 151


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--r1', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    votes = defaultdict(lambda: [Counter() for _ in range(BLOCK_FROM_R1)])
    depth = Counter()
    op = gzip.open if args.r1.endswith('.gz') else open
    with op(args.r1, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            block = line.strip()[:BLOCK_FROM_R1]
            if len(block) < 16 or not set(block[:16]) <= set('ACGT'):
                continue
            seed = int.from_bytes(dna_to_bytes(block[:16]), 'big')
            depth[seed] += 1
            row = votes[seed]
            for p, b in enumerate(block):
                if b in 'ACGT':
                    row[p][b] += 1
    with open(args.out, 'w') as f:
        for seed, pos_counters in votes.items():
            cons = ''.join(c.most_common(1)[0][0] if c else 'N' for c in pos_counters)
            minpos = min((sum(c.values()) for c in pos_counters if c), default=0)
            f.write(f"{seed}\t{cons}\t{depth[seed]}\t{minpos}\n")
    print(f"seeds with consensus: {len(votes)}", file=sys.stderr)


if __name__ == '__main__':
    main()
