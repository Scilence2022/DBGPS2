#!/usr/bin/env python3
"""Extend the block-12 parameter sweep without changing historical outputs."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import run_block_decode as runner

OUT = runner.DATA / 'results/block12_small_k'
READS = runner.OUT / 'reads/broad/block12.fq.gz'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['initial', 'fine', 'fine14'], default='initial')
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    historical = json.loads((runner.DATA / 'results/inner_rs/file3_candidate_path_bounds.json').read_text())
    previous = {(r['k'], r['cutoff'], float(r['ratio'])) for r in historical}
    grids = {
        16: [(c, r) for r in [1.5, 1.8, 2.1, 2.55, 3.0, 3.5, 5.0, 10.0] for c in [1, 2, 3, 5, 9]],
        14: [(c, r) for r in [1.2, 1.5, 1.8, 2.1, 2.55, 3.5] for c in [1, 3, 5, 9, 13, 19, 29]],
        15: [(c, r) for r in [1.8, 2.1, 2.55, 3.0, 3.5, 5.0] for c in [2, 4, 5, 7]],
        17: [(c, r) for r in [2.55, 3.0, 3.5, 4.0, 5.0, 10.0] for c in [1, 2, 3]],
        13: [(c, r) for r in [1.5, 2.1, 3.5] for c in [13, 19, 29]],
    }
    if args.stage == 'fine':
        # Follow up on the initial k16/c1/r3.5 and k15/c2/r2.1 maxima.
        grids = {
            16: [(1, r) for r in [2.8, 3.1, 3.2, 3.3, 3.4, 3.6, 3.75, 4.0, 4.25, 4.5]],
            15: [(c, r) for r in [1.9, 2.0, 2.2, 2.3, 2.4] for c in [1, 2]],
        }
    elif args.stage == 'fine14':
        # Initial k14 prefers c9; check intermediate cutoffs around its branching transition.
        grids = {14: [(c, r) for r in [1.8, 2.1, 2.55, 3.0, 3.5] for c in [6, 7, 8, 10, 11]]}
    grids = {k: [(c, r) for c, r in grid if (k, c, r) not in previous] for k, grid in grids.items()}
    manifest = dict(input=str(READS), input_reads=1681943, input_sha256=hashlib.sha256(READS.read_bytes()).hexdigest(),
                    global_index_domain=196596, block=12, block_size=16383, read_limit=200,
                    live_path_limit=1050, reference_assisted_read_assignment=True,
                    reference_use='Design assigns reads; graph search and inner-RS decisions are reference-blind; truth scores finalized outputs.',
                    plans=[dict(k=k, grid=grid) for k, grid in grids.items()],
                    new_configurations=sum(map(len, grids.values())),
                    baseline_configurations=179,
                    binary_sha256=hashlib.sha256((runner.DATA / 'bin/graph_block_rs').read_bytes()).hexdigest())
    initial = OUT / 'plan.initial.json'
    if args.stage != 'initial' and not initial.exists():
        initial.write_text((OUT / 'plan.json').read_text())
    (OUT / f'plan.{args.stage}.json').write_text(json.dumps(manifest, indent=2))
    merged = dict(manifest); combined = {}
    for path in sorted(OUT.glob('plan.*.json')):
        for row in json.loads(path.read_text())['plans']:
            combined.setdefault(row['k'], []).extend(row['grid'])
    merged['plans'] = [dict(k=k, grid=grid) for k, grid in sorted(combined.items())]
    merged['new_configurations'] = sum(len(grid) for grid in combined.values())
    (OUT / 'plan.json').write_text(json.dumps(merged, indent=2))
    runner.OUT = OUT
    with concurrent.futures.ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda pair: runner.job('broad' if args.stage == 'initial' else args.stage, 12, pair[0], pair[1], READS), grids.items()))
    print('SMALL-K SWEEP COMPLETE', flush=True)


if __name__ == '__main__':
    main()
