#!/usr/bin/env python3
"""Reference-scored parameter search; sequence generation never reads references.

All graph construction uses 100% of the raw reads. A fixed hash subset of observed
seeds screens parameters cheaply; selected settings are confirmed on ALL seeds.
References select/report parameters, never repair sequences or identify which
candidate to emit for an address. This is dataset-specific tuning, not a held-out
generalization estimate or a proof of a global mathematical optimum.
"""
from snapshot_provenance import snapshot_commit, snapshot_modified
from collections import defaultdict
import concurrent.futures as cf
import csv
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OPT = ROOT / "optimization"
REPO = ROOT.parents[1]
NOHP = re.compile("AA|CC|GG|TT")
TRUTH = set()
BY_PREFIX = {}
with gzip.open(ROOT / "data/id20.refs.txt.gz", "rt") as handle:
    for line in handle:
        body = line.strip()[3:-3][20:130]
        TRUTH.add(body)
        BY_PREFIX[body[:14]] = body
assert len(TRUTH) == len(BY_PREFIX) == 607150


def score_file(path, detail=False):
    """Reference-free conflict rejection, followed by reference-based scoring."""
    bodies = {"F": {}, "R": {}, "both": {}}
    source_rows = 0
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            body = row["payload"]
            assert len(body) == 110
            if NOHP.search(body):
                continue
            source_rows += 1
            prefix = body[:14]
            for scope in [row["direction"], "both"]:
                old = bodies[scope].get(prefix)
                if old is None:
                    bodies[scope][prefix] = body
                elif old != body:
                    # Retain a permanent conflict marker; never choose by truth.
                    bodies[scope][prefix] = "CONFLICT"
    result = {"source_rows_after_nohp": source_rows}
    for scope, candidates in bodies.items():
        accepted = [body for body in candidates.values() if body != "CONFLICT"]
        exact = sum(body in TRUTH for body in accepted)
        known_wrong = sum(body not in TRUTH and body[:14] in BY_PREFIX for body in accepted)
        result[scope] = {"exact": exact, "wrong": len(accepted)-exact, "wrong_known_address": known_wrong,
                         "unknown_address": len(accepted)-exact-known_wrong,
                         "returned": len(accepted), "conflict_addresses": sum(body == "CONFLICT" for body in candidates.values()),
                         "precision": exact/len(accepted) if accepted else 0,
                         "exact_rate": exact/607150}
    return (result, bodies) if detail else result


def tag_for(k, cutoff, ratio, skip):
    return f"k{k}_c{cutoff}_r{ratio:.3f}_s{skip}"


class Graph:
    def __init__(self, k, phase):
        self.k, self.phase = k, phase
        self.log = (OPT / "results" / f"{phase}.k{k}.core.log").open("w")
        self.started = time.monotonic()
        binary = "search_core_inner" if phase.startswith("inner") else "search_core"
        self.process = subprocess.Popen([str(OPT / "tmp" / binary), str(ROOT / "data/full.fa"), str(k), "4"],
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.log, text=True, bufsize=1)
        if self.process.stdout.readline().strip() != "READY":
            raise RuntimeError(f"Graph k={k} failed")
        print(f"{phase}: k={k} graph ready ({time.monotonic()-self.started:.1f}s)", flush=True)

    def run(self, cutoff, ratio, skip, seeds):
        tag = tag_for(self.k, cutoff, ratio, skip)
        out = OPT / ("tmp" if self.phase.endswith("screen") else "results") / f"{self.phase}.{tag}.tsv"
        started = time.monotonic()
        self.process.stdin.write(f"{tag} {cutoff} {ratio:.8f} {skip} {seeds} {out}\n")
        self.process.stdin.flush()
        response = self.process.stdout.readline().split()
        if not response or response[0] != "DONE" or response[1] != tag:
            raise RuntimeError(f"Search failed: {tag}, {response}")
        result = {"tag": tag, "k": self.k, "cutoff": cutoff, "ratio": ratio, "skip": skip,
                  "queries": int(response[2]), "raw_outputs": int(response[3]), "search_cpu_seconds": float(response[4]),
                  "wall_seconds": time.monotonic()-started, "output": str(out), "phase": self.phase}
        result.update(score_file(out))
        print(f'{self.phase}: {tag}, F={result["F"]["exact"]}, R={result["R"]["exact"]}, '
              f'both={result["both"]["exact"]}, errors={result["both"]["wrong"]}, {result["wall_seconds"]:.1f}s', flush=True)
        with (OPT / "results" / f"{self.phase}.k{self.k}.jsonl").open("a") as handle:
            handle.write(json.dumps(result)+"\n")
        return result

    def close(self):
        self.process.stdin.write("QUIT\n"); self.process.stdin.flush()
        self.process.stdin.close()
        assert self.process.wait() == 0
        self.log.close()


def screen_k(k, phase="screen"):
    graph = Graph(k, phase)
    results = []
    done = set()
    def run(c, r, s):
        key = c, round(r, 3), s
        if key in done:
            return
        done.add(key)
        results.append(graph.run(*key, OPT / "tmp/screen.seeds.tsv"))
    for c, r in itertools.product([0, 1, 2, 3, 4, 6], [1.1, 1.25, 1.4, 1.6, 1.8, 2.0, 2.55, 3.0]):
        run(c, r, 0)
    best = sorted(results, key=lambda x: (x["both"]["exact"], -x["both"]["wrong"]), reverse=True)[:2]
    for row in best:
        for ratio in [row["ratio"]-.1, row["ratio"]-.05, row["ratio"]+.05, row["ratio"]+.1]:
            if ratio >= 1.05:
                run(row["cutoff"], ratio, 0)
        for skip in [4, 16]:
            run(row["cutoff"], row["ratio"], skip)
    graph.close()
    return results


def confirm_k(item):
    k, configs = item
    graph = Graph(k, "full")
    results = [graph.run(row["cutoff"], row["ratio"], row["skip"], OPT / "tmp/all.seeds.tsv") for row in configs]
    graph.close()
    return results


def main():
    started = time.monotonic()
    seed_meta = OPT / "results/seed_provenance.json"
    while not seed_meta.exists():
        time.sleep(2)
    cmd = ["cc", "-O3", "-Wall", "-I.", "-o", str(OPT / "tmp/search_core"),
           str(OPT / "scripts/search_core.c"), "dbgps_codec.c", "dbgps_redundancy.c", "dbgps_strand.c", "kthread.c", "-lz", "-lpthread"]
    subprocess.run(cmd, cwd=REPO, check=True)
    with cf.ThreadPoolExecutor(2) as pool:
        screened = [row for group in pool.map(screen_k, [31, 29, 27, 25, 23, 21]) for row in group]
    (OPT / "results/screen_summary.json").write_text(json.dumps(screened, indent=2))
    ranked = sorted(screened, key=lambda x: (x["both"]["exact"], -x["both"]["wrong"]), reverse=True)
    selected = {x["tag"]: x for x in ranked[:3]}
    precise = [x for x in ranked if x["both"]["precision"] >= .9999]
    for row in precise[:2]:
        selected[row["tag"]] = row
    # Include orientation-specific winners and the previous full-data baseline.
    for scope in ["F", "R"]:
        row = max(screened, key=lambda x: (x[scope]["exact"], -x[scope]["wrong"]))
        selected[row["tag"]] = row
    baseline = next(x for x in screened if x["k"] == 31 and x["cutoff"] == 1 and x["ratio"] == 2.55 and x["skip"] == 0)
    selected[baseline["tag"]] = baseline
    (OPT / "results/selected_configs.json").write_text(json.dumps(list(selected.values()), indent=2))
    groups = defaultdict(list)
    for row in selected.values():
        groups[row["k"]].append(row)
    print("Full-seed confirmations:", list(selected), flush=True)
    with cf.ThreadPoolExecutor(2) as pool:
        confirmed = [row for group in pool.map(confirm_k, sorted(groups.items(), reverse=True)) for row in group]
    (OPT / "results/full_summary.json").write_text(json.dumps(confirmed, indent=2))
    meta = {"read_fraction": 1.0, "raw_reads": 15719903, "raw_mean_coverage": 15719903/607150,
            "graph_source_commit": snapshot_commit(REPO),
            "production_source_modified": bool(snapshot_modified(REPO)),
            "source_sha256": hashlib.sha256((REPO/"DBGPS2.c").read_bytes()).hexdigest(),
            "screened_configs": len(screened), "full_confirmed_configs": len(confirmed), "elapsed_seconds": time.monotonic()-started,
            "references_used": "parameter scoring/selection and final exactness assessment only; never candidate generation/repair", "rs_validation": "not_performed"}
    (OPT / "results/provenance.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
