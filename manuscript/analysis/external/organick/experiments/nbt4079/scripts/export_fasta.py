#!/usr/bin/env python3
"""Export predicted 110 nt bodies after the reference-free homopolymer filter.

These files contain algorithm predictions, including their residual errors.
Reference sequences are not read, and no output is repaired using ground truth.
"""
import csv
import gzip
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
pattern = re.compile(r"AA|CC|GG|TT")
counts = {}
for sample in ["5x_raw", "10x_raw", "full"]:
    count = 0
    with (ROOT / "results" / (sample+".c1_r2.55.tsv")).open() as source, gzip.open(
            ROOT / "results" / (sample+".c1_r2.55.nohp.fa.gz"), "wt", compresslevel=1) as dest:
        for row in csv.DictReader(source, delimiter="\t"):
            body = row["payload"]
            assert len(body) == 110
            if pattern.search(body):
                continue
            dest.write(f'>seed_{row["seed"][-14:]} support={row["support"]}\n{body}\n')
            count += 1
    counts[sample] = count
summary = {row["run"]: row for row in json.loads((ROOT / "results/summary.json").read_text())}
for sample, count in counts.items():
    assert count == summary[sample+".c1_r2.55.nohp"]["returned_outputs"]
(ROOT / "results/fasta_exports.json").write_text(json.dumps(counts, indent=2))
print(json.dumps(counts))
