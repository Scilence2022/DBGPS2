#!/usr/bin/env python3
"""Revisit File 3 configurations for single-parameter maxima or parameter unions."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import subprocess
import time

ROOT=pathlib.Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/inner_rs'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parallel',type=int,default=3)
    ap.add_argument('--k',type=int,nargs='+',default=[19,20,21,18,17,15,25,31,13])
    ap.add_argument('--remaining',action='store_true',help='Run all missing original configurations, without the single-configuration bound pruning')
    args=ap.parse_args()
    bounds=json.loads((OUT/'file3_candidate_path_bounds.json').read_text())
    def run(k):
        def completed(r):
            stem=f"file3.k{k}.c{r['cutoff']}.r{r['ratio']:.2f}"
            # Metadata is written before the diagnostic, so also check its final row.
            diag=OUT/(stem+'.truth_candidates.tsv')
            if not (OUT/(stem+'.inner_stats.json')).exists() or not diag.exists():return False
            with diag.open('rb') as handle:
                handle.seek(max(0,diag.stat().st_size-256));lines=handle.read().splitlines()
            last=lines[-1].split(b'\t') if lines else []
            return len(last)==4 and last[0]==b'196595'
        grid=[r for r in bounds if r['k']==k and (not completed(r) if args.remaining else not r['dominated_by_strict_baseline'])]
        grid.sort(key=lambda r:(r['ratio'],r['cutoff']))
        if not grid:
            print('SKIP',k,'no missing configurations',flush=True);return
        suffix='_union_remaining' if args.remaining else ''
        print('START',k,'configurations',len(grid),flush=True)
        env=os.environ.copy();env['ANTK_LAYOUT']='internal';env['ANTK_BLOOM']='1'
        env['ANTK_GRID']=','.join(f"{r['ratio']:.2f}:{r['cutoff']}" for r in grid)
        cmd=[str(DATA/'bin/graph_inner_rs'),str(DATA/'processed/file3.fa.gz'),'196596',str(k),str(OUT/'file3'),str(DATA/'sources/File3_synthesized.fasta')]
        start=time.monotonic()
        with (OUT/f'file3_k{k}{suffix}.log').open('w') as log:
            p=subprocess.run(['/usr/bin/time','-l']+cmd,env=env,stdout=log,stderr=log)
        record=dict(k=k,command=cmd,grid=grid,exit_code=p.returncode,seconds=time.monotonic()-start)
        controls=[]
        if not p.returncode:
            for r in grid:
                name=f"file3.k{k}.c{r['cutoff']}.r{r['ratio']:.2f}"
                strict=OUT/(name+'.strands');collector=OUT/(name+'.strict.strands');old=DATA/'results'/(name+'.strands')
                assert strict.read_bytes()==old.read_bytes()==collector.read_bytes(),name
                controls.append(dict(name=name,byte_identical=True,sha256=hashlib.sha256(strict.read_bytes()).hexdigest()))
        record['strict_controls']=controls
        (OUT/f'file3_k{k}{suffix}.run.json').write_text(json.dumps(record,indent=2))
        print('FINISHED',k,'exit',p.returncode,'strict controls identical',len(controls),'seconds',round(record['seconds'],2),flush=True)
        if p.returncode:raise RuntimeError(record)
    with concurrent.futures.ThreadPoolExecutor(args.parallel) as pool:list(pool.map(run,args.k))

if __name__=='__main__':main()
