#!/usr/bin/env python3
"""Restore the synthesized wire layout from later deposited design codewords.

The later File 2/3 design files store (in reverse base order) an unrandomized
prefix index, payload, and freshly generated parity. The physical experiment
uses payload | randomized index | parity, as in the primary File 1 reference.
Only evaluation references are transcoded. Assembly receives numeric indices,
not reference payloads. All File 1 strands verify this transformation exactly.
"""
import pathlib
import json
import hashlib

ROOT = pathlib.Path(__file__).resolve().parents[2] / 'input/antkowiak2020'
ALPHABET = 'ACGT'
INDEX_XOR = 1071068

def mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 64:
            a ^= 91
    return result

def encode_word(value, bases):
    return ''.join(ALPHABET[(value >> (2*i)) & 3] for i in range(bases-1, -1, -1))

def value(s):
    n = 0
    for c in s:
        n = (n << 2) | ALPHABET.index(c)
    return n

INV_THREE = next(x for x in range(64) if mul(x, 3) == 1)

def synthesized(payload, index):
    prefix = payload + encode_word(index ^ INDEX_XOR, 12)
    s0 = s1 = 0
    power = mul(2, 2)
    for j in range(18):
        symbol = value(prefix[51-3*j:54-3*j])
        s0 ^= symbol
        s1 ^= mul(symbol, power)
        power = mul(power, 2)
    c1 = mul(s0 ^ s1, INV_THREE)
    c0 = s0 ^ c1
    return prefix + encode_word(c1, 3) + encode_word(c0, 3)

def read_fasta(path):
    return [s for s in path.read_text().splitlines() if s and not s.startswith('>')]

def main():
    primary = (ROOT/'sources/noisy_dna_data_storage/data/File1_ODNA.txt').read_text().splitlines()
    assert all(synthesized(s[:42], i) == s for i, s in enumerate(primary))
    manifest = {'primary_File1_reencoded_exact': len(primary), 'index_xor': INDEX_XOR,
                'inner_field_polynomial': 91, 'inner_generator_roots': [1, 2]}
    for label in [2, 3]:
        source = ROOT/f'sources/File{label}_design.fasta'
        original = read_fasta(source)
        target = source.with_name(f'File{label}_synthesized.fasta')
        with target.open('w') as out:
            for i, s in enumerate(original):
                normalized = s[::-1]
                assert value(normalized[:12]) == i
                out.write(f'>{i:06d}\n{synthesized(normalized[12:54], i)}\n')
        manifest[f'File{label}'] = dict(original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                                      evaluation_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                                      sequences=len(original), transform='reverse bases; preserve 42 nt payload; reserialize randomized index and inner RS parity')
    (ROOT/'results/reference_transcoding.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
