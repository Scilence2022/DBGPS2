#!/usr/bin/env python3
"""Build adapter, run three nested read sets, and score exact reconstruction."""
import concurrent.futures as cf
import hashlib
import json
from pathlib import Path
import platform
import resource
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]


def run(sample):
    started = time.monotonic()
    cmd = [str(ROOT / "tmp/dbgps2_adapter"), str(ROOT / "data" / (sample+".fa")),
           str(ROOT / "data" / (sample+".seeds.tsv")), str(ROOT / "results" / sample), "4"]
    with (ROOT / "results" / (sample+".stdout")).open("w") as out, (ROOT / "results" / (sample+".log")).open("w") as log:
        subprocess.run(cmd, stdout=out, stderr=log, check=True)
    result = {"sample": sample, "command": cmd, "wall_seconds": time.monotonic()-started,
              "completed_configs": ["c0_r2.55", "c1_r2.55"]}
    (ROOT / "results" / (sample+".timing.json")).write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    cmd = ["cc", "-O3", "-Wall", "-I.", "-o", str(ROOT / "tmp/dbgps2_adapter"),
           str(ROOT / "scripts/dbgps2_adapter.c"), "dbgps_codec.c", "dbgps_redundancy.c",
           "dbgps_strand.c", "kthread.c", "-lz", "-lpthread"]
    subprocess.run(cmd, cwd=REPO, check=True)
    source_files = ["DBGPS2.c", "dbgps_codec.c", "dbgps_redundancy.c", "dbgps_strand.c", "kthread.c"]
    modified = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", *source_files], cwd=REPO).returncode
    metadata = {"git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
                "compile_command": cmd, "platform": platform.platform(),
                "compiler": subprocess.check_output(["cc", "--version"], text=True).splitlines()[0],
                "sources_modified": bool(modified), "adapter": "scripts/dbgps2_adapter.c",
                "source_sha256": {name: hashlib.sha256((REPO/name).read_bytes()).hexdigest() for name in source_files},
                "graph_k": 31, "index_bases_internal": 16, "path_bases_internal": 96,
                "output": "2 observed forward-primer bases plus 110 nt encoded region; strip first 2 bases",
                "seed": "last 17 forward-primer bases + first 14 encoded bases, discovered only from reads",
                "redundancy": "none", "threads_per_graph": 4}
    metadata["completed_configs"] = ["c0_r2.55", "c1_r2.55"]
    (ROOT / "results/provenance.json").write_text(json.dumps(metadata, indent=2))
    # Independent processes avoid DBGPS2's global path-search scratch buffers.
    with cf.ThreadPoolExecutor(3) as pool:
        list(pool.map(run, ["5x_raw", "10x_raw", "full"]))
    subprocess.run(["python3", str(ROOT / "scripts/evaluate.py")], check=True)
