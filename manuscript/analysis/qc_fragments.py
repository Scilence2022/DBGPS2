#!/usr/bin/env python3
"""Fragment-length distribution from paired-end alignments.

For each read pair whose mates both align (primary) to the same designed strand,
the fragment length is the reference span covered by the two mates:
max(target_end) - min(target_start). Unlike the adapter position (which cannot
exceed the 150-nt read length) this resolves fragments up to the full 200-nt
strand. Pairs with only one aligned mate, or none, are counted separately.

usage: qc_fragments.py results/qc sample [sample ...]
writes results/qc/<sample>.frag.tsv (length, pairs) and prints a summary line.
"""
import sys, os, json
import numpy as np

qc = sys.argv[1]
MAXL = 200


def load(path):
    d = {}
    with open(path) as f:
        for line in f:
            t = line.split("\t", 13)
            if "tp:A:P" not in line:
                continue
            if t[0] in d:
                continue
            d[t[0]] = (t[5], int(t[7]), int(t[8]))
    return d


for s in sys.argv[2:]:
    a, b = load(f"{qc}/{s}_1.paf"), load(f"{qc}/{s}_2.paf")
    n_pairs = json.load(open(f"{qc}/{s}_1.json"))["reads"]
    h = np.zeros(MAXL + 1, dtype=np.int64)
    both = single = diff = 0
    for k, (t1, s1, e1) in a.items():
        if k in b:
            t2, s2, e2 = b[k]
            if t1 != t2:
                diff += 1
                continue
            h[min(max(e1, e2) - min(s1, s2), MAXL)] += 1
            both += 1
        else:
            single += 1
    single += sum(1 for k in b if k not in a)
    none = n_pairs - both - single - diff
    np.savetxt(f"{qc}/{s}.frag.tsv", np.c_[np.arange(MAXL + 1), h], fmt="%d", delimiter="\t",
               header="length\tpairs", comments="")
    cum = np.cumsum(h) / h.sum()
    q = lambda p: int(np.searchsorted(cum, p))
    print(json.dumps(dict(sample=s, pairs=n_pairs, both=both / n_pairs, one_mate=single / n_pairs,
                          different_strands=diff / n_pairs, none=none / n_pairs,
                          median=q(0.5), q25=q(0.25), q75=q(0.75),
                          frac_ge190=float(h[190:].sum() / h.sum()),
                          frac_lt100=float(h[:100].sum() / h.sum()))))
