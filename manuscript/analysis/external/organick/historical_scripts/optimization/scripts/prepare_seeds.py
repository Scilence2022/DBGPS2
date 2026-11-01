#!/usr/bin/env python3
"""Discover reverse seeds and sample query seeds without reading references."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]
OPT = ROOT / "optimization"
F = b"AGTGCAACAAGTCAATCCGT"
R = b"AATTGAATGCTTGCTTGCCG"
TRANS = bytes.maketrans(b"ACGTN", b"TGCAN")
RC_R = R.translate(TRANS)[::-1]
ANCHOR = RC_R[-17:]
reverse = Counter()
started = time.monotonic()
with (ROOT / "data/full.fa").open("rb") as reads:
    count = 0
    while header := reads.readline():
        seq = reads.readline().rstrip()
        assert header.startswith(b">") and seq
        for oriented in (seq, seq.translate(TRANS)[::-1]):
            pos = oriented.find(ANCHOR)
            if pos >= 0:
                prefix = oriented[pos+17:pos+31]
                if len(prefix) == 14 and all(c in b"ACGT" for c in prefix):
                    reverse[ANCHOR+prefix] += 1
        count += 1
        if count % 3000000 == 0:
            print(f"Reverse seed discovery: {count:,} reads, {time.monotonic()-started:.1f}s", flush=True)
sample_count = Counter(); all_count = Counter()
with (OPT / "tmp/all.seeds.tsv").open("wb") as full, (OPT / "tmp/screen.seeds.tsv").open("wb") as screen:
    for direction in [b"F", b"R"]:
        if direction == b"F":
            handle = (ROOT / "data/full.seeds.tsv").open("rb")
            items = ((line.split()[0], int(line.split()[1])) for line in handle)
        else:
            items = sorted(reverse.items())
        for seed, support in items:
            row = direction+b"\t"+seed+b"\t"+str(support).encode()+b"\n"
            full.write(row); all_count[direction.decode()] += 1
            # Fixed hash query sample, independent of reference labels/outcomes.
            if hashlib.blake2b(direction+seed, digest_size=8).digest()[0] < 2:
                screen.write(row); sample_count[direction.decode()] += 1
        if direction == b"F":
            handle.close()
metadata = {"all_queries": all_count, "screen_queries": sample_count, "sampling_probability": 1/128,
            "reference_sequences_used": False, "raw_records": count,
            "forward_anchor": F[-17:].decode(), "reverse_anchor": ANCHOR.decode(),
            "elapsed_seconds": time.monotonic()-started}
(OPT / "results/seed_provenance.json").write_text(json.dumps(metadata, indent=2))
print(json.dumps(metadata, indent=2), flush=True)
