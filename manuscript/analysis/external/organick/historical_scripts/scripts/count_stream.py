#!/usr/bin/env python3
"""Count and validate gzip while the range downloader is still running."""
from pathlib import Path
import json
import time
import zlib

root = Path(__file__).resolve().parents[1]
size = 1363110383
chunk_size = 2**20
decompressor = zlib.decompressobj(31)
lines = 0
raw_bytes = 0
started = time.monotonic()
for i in range((size + chunk_size - 1) // chunk_size):
    part = root / "tmp/id20.fastq.gz.chunks" / str(i)
    complete = root / "data/id20.fastq.gz"
    expected = min(chunk_size, size-i*chunk_size)
    while not (part.exists() and part.stat().st_size == expected) and not complete.exists():
        time.sleep(0.5)
    if part.exists():
        block = part.read_bytes()
    else:
        with complete.open("rb") as handle:
            handle.seek(i * chunk_size)
            block = handle.read(min(chunk_size, size - i * chunk_size))
    assert len(block) == min(chunk_size, size-i*chunk_size)
    while block:
        data = decompressor.decompress(block)
        lines += data.count(b"\n")
        raw_bytes += len(data)
        block = decompressor.unused_data
        if block:
            assert decompressor.eof
            decompressor = zlib.decompressobj(31)
    if i % 100 == 0:
        print(f"counted {lines//4:,} records through compressed chunk {i}", flush=True)
assert decompressor.eof
assert lines % 4 == 0
result = {"fastq_records": lines//4, "uncompressed_bytes": raw_bytes, "gzip_crc_validated": True}
(root / "results/fastq_count.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result), flush=True)
