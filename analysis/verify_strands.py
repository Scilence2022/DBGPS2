#!/usr/bin/env python3
"""Validate DBGPS2 *.strands reconstructions for DNA Fountain droplets.

Checks per line `seed\\tDNA152`:
- index bases decode back to the seed
- RS(38,36) validity (optionally corrected)
- biochemical screen (homopolymer <=3, GC 45-55%)
- exact agreement with the read-level position-wise consensus
"""
import argparse
import sys

from fountain_common import PRNG, ReedSolomon2, dna_to_bytes, screen_dna


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('strands')
    ap.add_argument('--consensus')
    ap.add_argument('--report_corrected', action='store_true')
    args = ap.parse_args()

    rs = ReedSolomon2()
    cons = {}
    if args.consensus:
        with open(args.consensus) as f:
            for line in f:
                seed, seq, _depth, _min = line.split('\t')
                cons[int(seed)] = seq

    n = index_mismatch = rs_ok = rs_corrected = screen_fail = 0
    exact = hamming_total = compared = 0
    mismatched = []
    with open(args.strands) as f:
        for line in f:
            seed_str, dna = line.strip().split('\t')
            seed = int(seed_str)
            n += 1
            if len(dna) != 152:
                index_mismatch += 1
                continue
            if int.from_bytes(dna_to_bytes(dna[:16]), 'big') != seed:
                index_mismatch += 1
                continue
            word = dna_to_bytes(dna)
            corrected = rs.check(word)
            if corrected is None:
                continue
            if corrected != list(word[:36]):
                if not args.report_corrected:
                    continue
                rs_corrected += 1
            rs_ok += 1
            if not screen_dna(dna[:136]):
                screen_fail += 1
            if seed in cons:
                compared += 1
                d = sum(a != b for a, b in zip(dna, cons[seed]))
                hamming_total += d
                if d == 0:
                    exact += 1
                elif len(mismatched) < 20:
                    mismatched.append((seed, d))
    print(f"reconstructed: {n}")
    print(f"index-consistent: {n - index_mismatch}")
    print(f"RS-valid{'(+corrected)' if args.report_corrected else ''}: {rs_ok} "
          f"(corrected: {rs_corrected})")
    print(f"screen failures among RS-valid: {screen_fail}")
    if compared:
        print(f"consensus-compared: {compared}, exact matches: {exact} "
              f"({100.0 * exact / compared:.2f}%), total mismatches: {hamming_total}")
        if mismatched:
            print("first mismatched seeds:", mismatched)
    return 0


if __name__ == '__main__':
    sys.exit(main())
