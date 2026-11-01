#!/usr/bin/env python3
"""Download author data at an immutable revision and verify Git LFS SHA256."""
import concurrent.futures as cf
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
REV = "1077c74a0e002be83ab832e69bfd560e94d4d202"
BASE = f"https://media.githubusercontent.com/media/uwmisl/data-nbt17/{REV}/"


def download(name, workers=32, chunk_size=2**20):
    pointer_url = f"https://raw.githubusercontent.com/uwmisl/data-nbt17/{REV}/{name}"
    with urllib.request.urlopen(pointer_url, timeout=60) as response:
        pointer = response.read().decode()
    (ROOT / "sources" / (name + ".lfs.txt")).write_text(pointer)
    fields = dict(line.split(" ", 1) for line in pointer.splitlines())
    expected = fields["oid"].split(":")[1]
    size = int(fields["size"])
    dest = ROOT / ("sources" if name == "MANIFEST" else "data") / name
    if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest() == expected:
        print(name, "already verified", flush=True)
        return
    parts = ROOT / "tmp" / (name + ".chunks")
    parts.mkdir(exist_ok=True)
    started = time.monotonic()

    def fetch(i):
        start, end = i * chunk_size, min((i + 1) * chunk_size, size) - 1
        part = parts / str(i)
        if part.exists() and part.stat().st_size == end - start + 1:
            return part.stat().st_size
        for attempt in range(5):
            try:
                request = urllib.request.Request(BASE + name, headers={"Range": f"bytes={start}-{end}"})
                with urllib.request.urlopen(request, timeout=90) as response:
                    if response.status != 206 and size > chunk_size:
                        raise ValueError(f"Range ignored: HTTP {response.status}")
                    content_range = response.headers.get("Content-Range", "")
                    if response.status == 206 and not content_range.startswith(f"bytes {start}-{end}/"):
                        raise ValueError(f"Wrong range: {content_range}")
                    data = response.read()
                if len(data) != end - start + 1:
                    raise ValueError(f"Length {len(data)} != {end-start+1}")
                part.write_bytes(data)
                return len(data)
            except Exception as exc:
                if attempt == 4:
                    raise
                time.sleep(1 + attempt)

    done = 0
    last = started
    with cf.ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(fetch, i) for i in range((size + chunk_size - 1) // chunk_size)]
        for future in cf.as_completed(futures):
            done += future.result()
            now = time.monotonic()
            if now - last > 20 or done == size:
                print(f"{name}: {done:,}/{size:,} bytes ({done/size:.1%}), {done/(now-started)/1e6:.2f} MB/s", flush=True)
                last = now
    sha = hashlib.sha256()
    with dest.with_suffix(dest.suffix + ".verified-part").open("wb") as output:
        for i in range(len(futures)):
            data = (parts / str(i)).read_bytes()
            sha.update(data)
            output.write(data)
    if sha.hexdigest() != expected:
        raise ValueError(f"SHA256 mismatch for {name}")
    os.replace(dest.with_suffix(dest.suffix + ".verified-part"), dest)
    metadata = {"file": name, "url": BASE + name, "bytes": size, "sha256": expected,
                "verified_against_git_lfs": True, "seconds": time.monotonic() - started}
    (ROOT / "sources" / (name + ".download.json")).write_text(json.dumps(metadata, indent=2))
    print(metadata, flush=True)
    for part in parts.iterdir():
        part.unlink()
    parts.rmdir()


if __name__ == "__main__":
    import sys
    for name in sys.argv[1:] or ["MANIFEST", "id20.refs.txt.gz", "id20.fastq.gz"]:
        download(name)
