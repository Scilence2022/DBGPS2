#!/usr/bin/env python3
"""Per-block table of photolithographic File 3 behind Fig. 4c-f.

Joins read assignment, k-mer coverage, pooled and routed outputs, the 179-setting
candidate union and the loss causes at the best pooled setting (k19/c1/r3.5), and
summarises correct-candidate coverage per tenth of each block in design order.
"""
import csv
import gzip
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, "external", "photolithographic")
OUT = os.path.join(HERE, "data", "source_data", "photolith_file3_blocks.csv")


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


report = open(os.path.join(EV, "reports", "File3.en.md")).read()
pooled_sec = report[report.index("## 6. Before partitioning"):report.index("## 7. Candidate-level")]
union_sec = report[report.index("## 7."):report.index("## 8.")]
pooled = {int(m[0]): m[1:] for m in re.findall(
    r"^\| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| (?:True|False) \|$", pooled_sec, flags=re.M)}
union = {int(b): u for b, u in re.findall(r"^\| (\d+) \| (\d+) \| \d+ \| (?:True|False) \|$", union_sec, flags=re.M)}
routed = {x["block"]: x for x in json.load(open(os.path.join(HERE, "data", "source_data", "figure_data.json")))["photo_selected_blocks"]}
reads = {b["block"]: b["broad_reads"] for b in json.load(open(os.path.join(EV, "evidence", "File3", "block_decode", "reads", "split_counts.json")))["blocks"]}
cover = {int(r["block"]): r for r in tsv(os.path.join(EV, "evidence", "File3", "analyzer_block12", "all_blocks.tsv")) if r["dataset"] == "processed.k19.L60"}
loss = {int(r["block"]): r for r in tsv(os.path.join(EV, "evidence", "File3", "analyzer_block12", "best_parameter_loss_breakdown.tsv"))}

fields = ["block", "assigned_reads", "median_kmer_coverage_k19", "missing_window_fraction", "strands_all_kmers_absent",
          "pooled_correct", "pooled_wrong", "pooled_erased", "pooled_budget", "routed_correct", "routed_wrong", "routed_erased",
          "routed_budget", "union_179_correct", "failed_candidates", "absent_kmer", "singleton_only", "ratio_rejected", "overflow"]
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(fields)
    for b in range(1, 13):
        c, l, r = cover[b], loss[b], routed[b]
        w.writerow([b, reads[b], round(float(c["median_strand_mean_coverage"]), 3), round(float(c["missing_window_fraction"]), 5),
                    c["all_kmers_zero"], *pooled[b], r["correct"], r["wrong"], r["erasures"], r["budget"], union[b],
                    l["failed_after_rs"], l["missing_kmer"], l["singleton_bottleneck"], l["rounded_ratio_rejected"], l["feasible_but_overflow"]])
print("wrote", OUT)

STRANDS = 16383
per_index = os.path.join(EV, "evidence", "File3", "inner_rs", "parameter_union", "per_index.tsv.gz")
counts = {}
for r in csv.DictReader(gzip.open(per_index, "rt"), delimiter="\t"):
    key = (int(r["block_1based"]), min(9, int(r["offset_0based"]) * 10 // STRANDS))
    n, ok = counts.get(key, (0, 0))
    counts[key] = (n + 1, ok + int(r["candidate_success"]))
position = OUT.replace("_blocks.csv", "_position.csv")
with open(position, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["block", "tenth", "strands", "correct_candidates", "percent"])
    for (b, d), (n, ok) in sorted(counts.items()):
        w.writerow([b, d + 1, n, ok, round(100 * ok / n, 2)])
print("wrote", position)
