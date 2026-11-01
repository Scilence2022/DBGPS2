#!/usr/bin/env python3
"""Uniform nested sampling of raw reads and reference-free discovery of seeds.

No reference FASTA is opened by this program. All raw reads, including N-bearing
and incomplete reads, go into the graph. Forward primer suffix matching is used
only to discover seed candidates and a simple modal-read baseline.
"""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import random
import time

ROOT = Path(__file__).resolve().parents[1]
FORWARD = b"AGTGCAACAAGTCAATCCGT"
ANCHOR = FORWARD[-17:]
COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")
BASES = set(b"ACGT")


def iter_reads(path):
    with gzip.open(path, "rb") as handle:
        while header := handle.readline():
            seq, plus, qual = handle.readline().rstrip(), handle.readline(), handle.readline().rstrip()
            if not header.startswith(b"@") or not plus.startswith(b"+") or len(seq) != len(qual):
                raise ValueError("Invalid FASTQ record")
            yield seq


def prepare(total):
    # 607150 is the published number of designed strands, not a read-derived label.
    targets = {"5x_raw": 5 * 607150, "10x_raw": 10 * 607150, "full": total}
    probabilities = {key: min(1, value / total) for key, value in targets.items()}
    outputs = {key: (ROOT / "data" / (key + ".fa")).open("wb") for key in targets}
    seeds = {key: Counter() for key in targets}
    modes = {key: {} for key in targets}
    stats = {key: Counter() for key in targets}
    rng = random.Random(20261017)
    started = time.monotonic()
    length_counts = Counter()
    for i, seq in enumerate(iter_reads(ROOT / "data/id20.fastq.gz"), 1):
        length_counts[len(seq)] += 1
        u = rng.random()
        selected = [key for key, probability in probabilities.items() if u < probability]
        observed = []
        for oriented in (seq, seq.translate(COMPLEMENT)[::-1]):
            start = oriented.find(ANCHOR)
            if start < 0:
                continue
            body_start = start + len(ANCHOR)
            prefix = oriented[body_start:body_start + 14]
            if len(prefix) != 14 or not set(prefix) <= BASES:
                continue
            seed = ANCHOR + prefix
            body = oriented[body_start:body_start + 110]
            # Include all 110 positions (no primer/inner-code/quality validation).
            complete = len(body) == 110 and set(body) <= BASES
            observed.append((seed, body if complete else None))
        for key in selected:
            outputs[key].write(b">" + str(i).encode() + b"\n" + seq + b"\n")
            stats[key]["raw_reads"] += 1
            stats[key]["N_bearing_reads"] += b"N" in seq
            stats[key]["anchored_reads"] += bool(observed)
            for seed, body in observed:
                seeds[key][seed] += 1
                if body is not None:
                    stats[key]["complete_anchored_reads"] += 1
                    modes[key].setdefault(seed, Counter())[body] += 1
        if i % 1000000 == 0:
            print(f"prepared {i:,} reads in {time.monotonic()-started:.1f}s", flush=True)
    if i != total:
        raise ValueError(f"Expected {total} records, observed {i}")
    for key in targets:
        outputs[key].close()
        with (ROOT / "data" / (key + ".seeds.tsv")).open("wb") as handle:
            for seed, count in sorted(seeds[key].items()):
                handle.write(seed + b"\t" + str(count).encode() + b"\n")
        with (ROOT / "results" / (key + ".modal.tsv")).open("wb") as handle:
            handle.write(b"seed\tsupport\tpayload\n")
            for seed, counter in sorted(modes[key].items()):
                common = counter.most_common(2)
                if len(common) == 1 or common[0][1] > common[1][1]:
                    body, count = common[0]
                    handle.write(seed + b"\t" + str(count).encode() + b"\t" + body + b"\n")
        stats[key]["distinct_seeds"] = len(seeds[key])
        stats[key]["mean_raw_coverage"] = stats[key]["raw_reads"] / 607150
    result = {"random_seed": 20261017, "lengths": dict(length_counts),
              "samples": stats, "elapsed_seconds": time.monotonic()-started,
              "reference_sequences_used_for_seeding": False}
    (ROOT / "results/preprocessing.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--total", type=int, required=True)
    args = parser.parse_args()
    prepare(args.total)
