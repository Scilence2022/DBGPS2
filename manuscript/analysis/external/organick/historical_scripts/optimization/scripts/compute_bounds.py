#!/usr/bin/env python3
"""Upper bounds from reference-path existence, explicitly separated from recovery."""
import concurrent.futures as cf
import json
import subprocess
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OPT=ROOT/"optimization"
REPO=ROOT.parents[1]
while not (OPT/"results/inner_full_summary.json").exists():
    time.sleep(2)
cmd=["cc","-O3","-Wall","-I.","-o",str(OPT/"tmp/graph_bound"),str(OPT/"scripts/graph_bound.c"),
     "dbgps_codec.c","dbgps_redundancy.c","dbgps_strand.c","kthread.c","-lz","-lpthread"]
subprocess.run(cmd,cwd=REPO,check=True)
def run(k):
    started=time.monotonic()
    p=subprocess.run([str(OPT/"tmp/graph_bound"),str(ROOT/"data/full.fa"),str(OPT/"tmp/all.seeds.tsv"),str(k),str(ROOT/"data/references.fa"),str(OPT/"tmp"/f"bound.k{k}.bin")],capture_output=True,text=True,check=True)
    row=json.loads(p.stdout);row['wall_seconds']=time.monotonic()-started
    (OPT/"results"/f"graph_bound.k{k}.json").write_text(json.dumps(row,indent=2))
    print(json.dumps(row),flush=True);return row
with cf.ThreadPoolExecutor(2) as pool:rows=list(pool.map(run,[31,29]))
(OPT/"results/graph_bounds.json").write_text(json.dumps(rows,indent=2))
maps=[(OPT/"tmp"/f"bound.k{k}.bin").read_bytes() for k in [31,29]]
assert len(maps[0])==len(maps[1])==607150
union=[a|b for a,b in zip(*maps)]
summary={"diagnostic":"oracle path-existence upper bound, not recovered strands", "k_union":[29,31],
         "bounds":[{"cutoff":c,"either_k_either_direction":sum((value>>c)&1 for value in union)} for c in range(7)]}
(OPT/"results/graph_bound_union.json").write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
