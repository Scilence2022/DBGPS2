#!/usr/bin/env python3
"""Optimize the paper-format candidate validator separately from the base run."""
import concurrent.futures as cf
from collections import defaultdict
import json
import subprocess
import time
from optimize import OPT, REPO, Graph, screen_k

started = time.monotonic()
while not (OPT / "results/full_summary.json").exists():
    time.sleep(2)
cmd = ["cc", "-O3", "-Wall", "-DUSE_INNER", "-I.", "-o", str(OPT / "tmp/search_core_inner"),
       str(OPT / "scripts/search_core.c"), "dbgps_codec.c", "dbgps_redundancy.c", "dbgps_strand.c", "kthread.c", "-lz", "-lpthread"]
subprocess.run(cmd, cwd=REPO, check=True)
with cf.ThreadPoolExecutor(2) as pool:
    screened = [row for group in pool.map(lambda k: screen_k(k, "inner_screen"), [31, 29]) for row in group]
(OPT / "results/inner_screen_summary.json").write_text(json.dumps(screened, indent=2))
ranked = sorted(screened, key=lambda row: (row["both"]["exact"], -row["both"]["wrong"]), reverse=True)
selected = {row["tag"]: row for row in ranked[:3]}
for scope in ["F", "R"]:
    row = max(screened, key=lambda x: (x[scope]["exact"], -x[scope]["wrong"]))
    selected[row["tag"]] = row
(OPT / "results/inner_selected_configs.json").write_text(json.dumps(list(selected.values()), indent=2))
groups = defaultdict(list)
for row in selected.values():
    groups[row["k"]].append(row)

def confirm(item):
    k, configs = item
    graph = Graph(k, "inner_full")
    results = [graph.run(row["cutoff"], row["ratio"], row["skip"], OPT / "tmp/all.seeds.tsv") for row in configs]
    graph.close()
    return results

with cf.ThreadPoolExecutor(2) as pool:
    confirmed = [row for group in pool.map(confirm, sorted(groups.items(), reverse=True)) for row in group]
(OPT / "results/inner_full_summary.json").write_text(json.dumps(confirmed, indent=2))
(OPT / "results/inner_provenance.json").write_text(json.dumps({"compile_command":cmd,
    "screened_configs":len(screened),"full_confirmed_configs":len(confirmed),
    "candidate_validator":"27 six-bit/ternary blocks plus zero terminal trit in AGCT rotation; no per-reference lookup",
    "production_graph_and_traversal_modified":False,"rs_validation":"not_performed",
    "elapsed_seconds_including_wait_for_base_run":time.monotonic()-started},indent=2))
print("Inner-code search complete", flush=True)
