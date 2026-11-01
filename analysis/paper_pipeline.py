#!/usr/bin/env python3
"""Paper-style read-level droplet pipeline (README of TeamErlich/dna-fountain):
merge each read pair (R1 + revcomp(R2), 150-nt overlap), keep exact-length
152-nt droplets, abundance-sort, and accept only RS-exact droplets
(receiver.py --max_hamming 0). Counts valid droplets and unique seeds."""
import argparse
import sys
from collections import Counter

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fountain_common import dna_to_bytes, ReedSolomon2, screen_dna

COMP = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(COMP)[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--r1', required=True)
    ap.add_argument('--r2', required=True)
    ap.add_argument('--out', required=True, help='valid droplet sequences, one per line')
    ap.add_argument('--max_overlap_mm', type=int, default=2,
                    help='max mismatches tolerated in the 150-nt overlap')
    args = ap.parse_args()
    rs = ReedSolomon2()

    merged = Counter()
    pairs = merged_pairs = 0
    with open(args.r1) as f1, open(args.r2) as f2:
        while True:
            h1 = f1.readline(); h2 = f2.readline()
            if not h1 or not h2:
                break
            r1 = f1.readline().strip()
            r2 = f2.readline().strip()
            pairs += 1
            if len(r1) < 151 or len(r2) < 151:
                continue
            if set(r1) > set('ACGT') or set(r2) > set('ACGT'):
                continue
            c2 = revcomp(r2)
            if set(r1[:16]) <= set('ACGT') and r1[1:151] == c2[:150]:
                # exact overlap merge; last base from R2
                merged[r1 + c2[150]] += 1
                merged_pairs += 1
    n_exact = sum(merged.values())
    valid = seen = screen_rej = 0
    seen_seeds = set()
    with open(args.out, 'w') as out:
        for seq, cnt in merged.most_common():
            word = dna_to_bytes(seq)
            dec = rs.check(word)
            if dec is None or dec != list(word[:36]):
                continue
            valid += 1
            seed = int.from_bytes(word[:4], 'big')
            if seed in seen_seeds:
                continue
            if not screen_dna(seq[:136]):
                screen_rej += 1
                continue
            seen_seeds.add(seed)
            out.write(seq + '\n')
            seen += 1
    print(f"pairs: {pairs}, exact-overlap merges: {n_exact} "
          f"({merged_pairs} distinct), RS-exact valid: {valid}, "
          f"screen-pass unique seeds: {seen} (screen rejects {screen_rej})",
          file=sys.stderr)


if __name__ == '__main__':
    main()
