#!/usr/bin/env python3
"""Peeling decoder for DNA Fountain droplets (glass.pyx port) with md5 report.

Input: file with one 152-nt droplet data block per line (ACGT), as produced by
DBGPS2 *.strands output after stripping the index column, or by consensus.
Only droplets that pass RS (optionally after correction) enter the decoder,
mirroring receiver.py --max_hamming 0 semantics by default.
"""
import argparse
import hashlib
import sys

sys.setrecursionlimit(10000000)

from fountain_common import PRNG, ReedSolomon2, dna_to_bytes, screen_dna


class Droplet:
    __slots__ = ('data', 'seed', 'num_chunks')

    def __init__(self, data, seed, num_chunks):
        self.data = data
        self.seed = seed
        self.num_chunks = set(num_chunks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile')
    ap.add_argument('--out', required=True)
    ap.add_argument('-n', '--chunk_num', type=int, default=67088)
    ap.add_argument('--header_size', type=int, default=4)
    ap.add_argument('--rs', type=int, default=2)
    ap.add_argument('--delta', type=float, default=0.001)
    ap.add_argument('--c_dist', type=float, default=0.025)
    ap.add_argument('--gc', type=float, default=0.05)
    ap.add_argument('-m', '--max_homopolymer', type=int, default=3)
    ap.add_argument('--no_screen', action='store_true',
                    help='skip the biochemical screen (for consensus-validated droplets)')
    ap.add_argument('--accept_corrected', action='store_true',
                    help='also accept droplets repaired by RS (default: exact only, max_hamming 0)')
    ap.add_argument('--max_line', type=int, default=None)
    args = ap.parse_args()

    rs = ReedSolomon2()
    prng = PRNG(K=args.chunk_num, delta=args.delta, c=args.c_dist)

    chunks = [None] * args.chunk_num
    done = set()
    chunk_to_droplets = [set() for _ in range(args.chunk_num)]
    seen_seeds = set()
    n_input = n_rs_valid = n_corrected = n_screen_reject = 0

    def update_entry(droplet):
        for cn in (droplet.num_chunks & done):
            seg = chunks[cn]
            droplet.data = bytes(droplet.data[i] ^ seg[i] for i in range(len(seg)))
            droplet.num_chunks.discard(cn)
            chunk_to_droplets[cn].discard(droplet)
        if len(droplet.num_chunks) == 1:
            lone = next(iter(droplet.num_chunks))
            chunks[lone] = droplet.data
            done.add(lone)
            droplet.num_chunks.clear()
            for other in list(chunk_to_droplets[lone]):
                update_entry(other)

    def add_droplet(droplet):
        for cn in droplet.num_chunks:
            chunk_to_droplets[cn].add(droplet)
        update_entry(droplet)

    with open(args.infile) as f:
        for line in f:
            if args.max_line and n_input >= args.max_line:
                break
            dna = line.strip().split('\t')[-1]
            if not dna or dna.startswith('>'):
                continue
            if len(dna) != 152:
                continue
            n_input += 1
            word = dna_to_bytes(dna)
            corrected = rs.check(word)
            if corrected is None:
                continue
            if corrected != list(word[:36]):
                if not args.accept_corrected:
                    continue
                n_corrected += 1
            n_rs_valid += 1
            seed = int.from_bytes(word[:4], 'big')
            payload = word[4:36]
            if seed in seen_seeds:
                continue
            seen_seeds.add(seed)
            _, _, nums = prng.get_src_blocks(seed)
            if not args.no_screen and not screen_dna(dna[:136], args.max_homopolymer, args.gc):
                n_screen_reject += 1
                continue
            add_droplet(Droplet(payload, seed, nums))

    sys.stderr.write(
        f"input droplets: {n_input}, RS-valid: {n_rs_valid} (corrected: {n_corrected}), "
        f"screen rejects: {n_screen_reject}, unique seeds: {len(seen_seeds)}, "
        f"chunks solved: {len(done)}/{args.chunk_num}\n")
    if len(done) < args.chunk_num:
        sys.stderr.write(f"DECODING INCOMPLETE: missing {args.chunk_num - len(done)} chunks\n")
    blob = b''.join(c if c is not None else b'\x00' * 32 for c in chunks)
    with open(args.out, 'wb') as f:
        f.write(blob)
    print(f"md5: {hashlib.md5(blob).hexdigest()}  size: {len(blob)}")
    print("expected md5: 8651e90d3a013178b816b63fdbb94b9b (2146816 bytes)")


if __name__ == '__main__':
    main()
