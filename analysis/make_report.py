#!/usr/bin/env python3
"""Aggregate DBGPS2 run metrics into a markdown comparison against the
DNA Fountain paper's reported numbers."""
import argparse
import hashlib
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from fountain_common import ReedSolomon2, dna_to_bytes


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_log(path):
    txt = open(path).read()
    def grab(pat, cast=float):
        m = re.search(pat, txt)
        return cast(m.group(1)) if m else None
    return {
        'kmer_time': grab(r'K-mer counting time: ([\d.]+)'),
        'decode_time': grab(r'Strand decoding time: ([\d.]+)'),
        'strands_decoded': grab(r'Strands decoded: (\d+)', int),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+', help='tag:dir pairs like onetile:analysis/runs/onetile')
    args = ap.parse_args()
    rows = []
    for item in args.runs:
        tag, d = item.split(':', 1)
        log = parse_log(f'{d}/dbgps.log')
        n_rs = n_rec = n_cons = n_exact = 0
        rs = ReedSolomon2()
        cons = {}
        try:
            for line in open(f'{d}/consensus.tsv'):
                seed, seq, *_ = line.split('\t')
                cons[int(seed)] = seq
        except FileNotFoundError:
            pass
        try:
            strands = open(f'{d}/dbgps.strands')
        except FileNotFoundError:
            strands = []
        for line in strands:
            seed, dna = line.strip().split('\t')
            n_rec += 1
            word = dna_to_bytes(dna)
            if rs.is_codeword(word):
                n_rs += 1
            if int(seed) in cons:
                n_cons += 1
                n_exact += all(a == b for a, b in zip(dna, cons[int(seed)]))
        digest = md5(f'{d}/decoded.bin') if __import__('os').path.exists(f'{d}/decoded.bin') else '-'
        rows.append((tag, n_rec, n_rs, log.get('decode_time'), digest))
    print('| run | DBGPS2 strands | RS-valid | decode time (s) | decoded md5 |')
    print('|---|---|---|---|---|')
    for tag, n_rec, n_rs, t, digest in rows:
        print(f'| {tag} | {n_rec} | {n_rs} | {t if t is not None else "-"} | {digest} |')
    print()
    print('Expected md5 (original file): 8651e90d3a013178b816b63fdbb94b9b')


if __name__ == '__main__':
    main()
