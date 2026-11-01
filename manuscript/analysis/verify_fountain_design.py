#!/usr/bin/env python3
"""Score DNA Fountain outputs against designs regenerated from the verified archive.

Re-runs the encoder (LFSR, robust soliton sampler, XOR, RS(38,36), biochemical
screen) on the MD5-verified 2,146,816-byte archive until 72,000 droplets pass,
then scores each output list by exact identity. Validity is tested as zero RS
syndrome with tables built here, independently of reedsolo.

Usage: verify_fountain_design.py FOUNTAIN_REPO OUT.csv [RUNS_DIR]
FOUNTAIN_REPO is a checkout of the dnafountain-adapter branch of DBGPS2 with
analysis/runs/{fullint,tile750k}; RUNS_DIR optionally holds benchmark
<case>.strands.txt files, one strand per line.
"""
import bisect
import csv
import glob
import hashlib
import os
import random
import sys

REPO, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(REPO, "analysis"))
from fountain_common import PRNG, bytes_to_dna, dna_to_bytes, py2_random_sample, screen_dna
from reedsolo import RSCodec

K, CHUNK, DESIGNS, MD5 = 67088, 32, 72000, "8651e90d3a013178b816b63fdbb94b9b"
LFSR_POLY, LFSR_START = 0b100000000000000000000000011000101, 0b001010101
RUNS = os.path.join(REPO, "analysis", "runs")

EXP, LOG = [0] * 512, [0] * 256
x = 1
for i in range(255):
    EXP[i], LOG[x] = x, i
    x <<= 1
    if x & 0x100:
        x ^= 0x11D
EXP[255:510] = EXP[:255]


def syndrome(word, root):
    v = 0
    for b in word:
        v = (EXP[LOG[v] + root] if v else 0) ^ b
    return v


def zero_syndrome(dna):
    word = dna_to_bytes(dna)
    return syndrome(word, 0) == 0 and syndrome(word, 1) == 0


def regenerate_designs(archive):
    blob = open(archive, "rb").read()
    assert hashlib.md5(blob).hexdigest() == MD5 and len(blob) == K * CHUNK
    chunks = [int.from_bytes(blob[i * CHUNK:(i + 1) * CHUNK], "big") for i in range(K)]
    cdf = PRNG(K=K, delta=0.001, c=0.025).cdf
    rs = RSCodec(2)
    nbits = LFSR_POLY.bit_length() - 1
    state, tried, designs = LFSR_START, 0, {}
    while len(designs) < DESIGNS:
        state <<= 1
        if (state >> nbits) & 1:
            state ^= LFSR_POLY
        tried += 1
        rnd = random.Random()
        rnd.seed(state)
        degree = min(bisect.bisect_right(cdf, rnd.random()) + 1, K)
        payload = 0
        for c in py2_random_sample(rnd, K, degree):
            payload ^= chunks[c]
        dna = bytes_to_dna(bytes(rs.encode(state.to_bytes(4, "big") + payload.to_bytes(CHUNK, "big"))))
        if screen_dna(dna):
            designs[state] = dna
    return designs, tried


def read_strands(path):
    return [line.split()[-1] for line in open(path) if line.strip()]


designs, tried = regenerate_designs(os.path.join(RUNS, "fullint", "decoded_do.bin"))
design_set = set(designs.values())
assert all(zero_syndrome(d) for d in design_set)
print(f"regenerated {len(designs):,} designs from {tried:,} LFSR states")

rows = []


def score(depth, method, path):
    strands = read_strands(path)
    valid = {d for d in strands if zero_syndrome(d)}
    rows.append({"depth": depth, "method": method, "reported": len(strands), "zero_syndrome": len(valid),
                 "exact_designs": len(valid & design_set), "not_design": len(valid - design_set)})
    return valid


graph = score("full", "DBGPS2 raw reads", os.path.join(RUNS, "fullint", "dbgps_do.strands"))
merge = score("full", "read merging", os.path.join(RUNS, "fullint", "paper_droplets.txt"))
graph1 = score("one_tile", "DBGPS2 raw reads", os.path.join(RUNS, "tile750k", "dbgps_fixed.strands"))
merge1 = score("one_tile", "read merging", os.path.join(RUNS, "tile750k", "paper_droplets.txt"))
if len(sys.argv) > 3:
    for f in sorted(glob.glob(os.path.join(sys.argv[3], "*.strands.txt"))):
        score("full", "benchmark " + os.path.basename(f).split(".")[0], f)

seeds = [int(s.split()[0]) for s in open(os.path.join(RUNS, "fullint", "seedlist.txt"))]
missed = set(designs) - {int.from_bytes(dna_to_bytes(d[:16]), "big") for d in graph}
summary = {"lfsr_states": tried, "queried_seeds": len(seeds), "queried_seeds_outside_design": sum(s not in designs for s in seeds),
           "designs_not_returned": len(missed), "designs_not_returned_unqueried": sum(s not in set(seeds) for s in missed),
           "shared_full": len(graph & merge), "graph_only_full": len(graph - merge), "merge_only_full": len(merge - graph),
           "shared_one_tile": len(graph1 & merge1), "graph_only_one_tile": len(graph1 - merge1),
           "merge_only_one_tile": len(merge1 - graph1), "union_one_tile": len(graph1 | merge1)}
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
with open(OUT.replace(".csv", "_summary.csv"), "w") as f:
    f.write("metric,value\n" + "".join(f"{k},{v}\n" for k, v in summary.items()))
for r in rows:
    print(r)
print(summary)
