#!/usr/bin/env python3
"""Alignment-based error profiling of subsampled reads.

Reads a minimap2 PAF (with --cs=long) for one read file plus the FASTQ, and
summarises per-sample error rates, substitution spectrum (oriented to the
designed strand), positional error profile, reads per strand, adapter
read-through and base quality.

usage: qc_align.py sample read_file.fq paf out_prefix
"""
import sys, re, json, collections
import numpy as np

sample, fq, paf, out = sys.argv[1:5]
STRAND_LEN = 200
ADAPTER = "AGATCGGAAGAGC"
COMP = str.maketrans("ACGTN", "TGCAN")

# ---- FASTQ: read count, quality, adapter read-through, poly-G ----
n_reads = 0; n_adapter = 0; n_polyg = 0; q_sum = 0; q_n = 0; q30 = 0
adapter_pos = collections.Counter()
with open(fq) as f:
    for i, line in enumerate(f):
        m = i % 4
        if m == 1:
            s = line.rstrip()
            n_reads += 1
            p = s.find(ADAPTER)
            if p >= 0:
                n_adapter += 1
                adapter_pos[p] += 1
            if "GGGGGGGGGGGGGGGGGGGG" in s:
                n_polyg += 1
        elif m == 3 and n_reads % 10 == 0:  # quality on 10% of reads
            q = np.frombuffer(line.rstrip().encode(), dtype=np.uint8) - 33
            q_sum += int(q.sum()); q_n += q.size; q30 += int((q >= 30).sum())

# ---- PAF ----
subs = collections.Counter()      # (ref_base, read_base) on designed strand
pos_mm = np.zeros(STRAND_LEN, dtype=np.int64)
pos_ins = np.zeros(STRAND_LEN, dtype=np.int64)
pos_del = np.zeros(STRAND_LEN, dtype=np.int64)
pos_cov = np.zeros(STRAND_LEN, dtype=np.int64)
start_hist = np.zeros(STRAND_LEN + 1, dtype=np.int64)
end_hist = np.zeros(STRAND_LEN + 1, dtype=np.int64)
per_strand = collections.Counter()
n_mapped = 0; aligned_bases = 0; n_mm = 0; n_ins = 0; n_del = 0; n_ins_ev = 0; n_del_ev = 0
n_perfect = 0; n_fullspan = 0
edit_hist = collections.Counter()
cs_re = re.compile(r"(=[ACGTN]+|\*[acgtn][acgtn]|\+[acgtn]+|-[acgtn]+)")
seen = set()
with open(paf) as f:
    for line in f:
        t = line.rstrip("\n").split("\t")
        qname = t[0]
        tags = t[12:]
        if not any(x == "tp:A:P" for x in tags):
            continue
        if qname in seen:
            continue
        seen.add(qname)
        n_mapped += 1
        strand = t[4]; tname = t[5]; ts = int(t[7]); te = int(t[8])
        per_strand[tname] += 1
        cs = next(x[5:] for x in tags if x.startswith("cs:Z:"))
        qlen = int(t[1]); qs = int(t[2]); qe = int(t[3])
        start_hist[ts] += 1; end_hist[te] += 1
        ed = 0
        rpos = ts
        for op in cs_re.findall(cs):
            c = op[0]
            if c == "=":
                L = len(op) - 1
                pos_cov[rpos:rpos + L] += 1
                aligned_bases += L
                rpos += L
            elif c == "*":
                rb, qb = op[1].upper(), op[2].upper()
                # cs is reported in target (designed strand) orientation; for
                # reverse-strand hits the read base is already complemented.
                subs[(rb, qb)] += 1
                pos_mm[rpos] += 1; pos_cov[rpos] += 1
                n_mm += 1; ed += 1; aligned_bases += 1
                rpos += 1
            elif c == "+":
                L = len(op) - 1
                n_ins += L; n_ins_ev += 1; ed += L
                if rpos < STRAND_LEN: pos_ins[rpos] += 1
            elif c == "-":
                L = len(op) - 1
                n_del += L; n_del_ev += 1; ed += L
                pos_del[rpos:rpos + L] += 1
                rpos += L
        edit_hist[min(ed, 20)] += 1
        if ed == 0:
            n_perfect += 1
        if ts <= 2 and te >= STRAND_LEN - 2:
            n_fullspan += 1

res = dict(sample=sample, file=fq, reads=n_reads, mapped=n_mapped,
           mapped_frac=n_mapped / max(n_reads, 1),
           adapter_frac=n_adapter / max(n_reads, 1), polyg_frac=n_polyg / max(n_reads, 1),
           mean_q=q_sum / max(q_n, 1), q30_frac=q30 / max(q_n, 1),
           aligned_bases=aligned_bases, mismatches=n_mm, ins_bases=n_ins, del_bases=n_del,
           ins_events=n_ins_ev, del_events=n_del_ev,
           sub_rate=n_mm / max(aligned_bases, 1), ins_rate=n_ins / max(aligned_bases, 1),
           del_rate=n_del / max(aligned_bases, 1),
           perfect_frac=n_perfect / max(n_mapped, 1),
           strands_hit=len(per_strand),
           subs={f"{a}>{b}": v for (a, b), v in subs.items()},
           edit_hist={int(k): v for k, v in edit_hist.items()},
           adapter_pos={int(k): v for k, v in adapter_pos.items()})
json.dump(res, open(out + ".json", "w"), indent=1)
np.savetxt(out + ".pos.tsv", np.c_[np.arange(STRAND_LEN), pos_cov, pos_mm, pos_ins, pos_del,
                                    start_hist[:STRAND_LEN], end_hist[:STRAND_LEN]],
           fmt="%d", delimiter="\t", header="pos\tcov\tmm\tins\tdel\tstart\tend", comments="")
with open(out + ".per_strand.tsv", "w") as f:
    for k, v in per_strand.items():
        f.write(f"{k}\t{v}\n")
print(json.dumps({k: res[k] for k in ("sample", "reads", "mapped_frac", "sub_rate", "ins_rate", "del_rate", "perfect_frac", "adapter_frac")}))
