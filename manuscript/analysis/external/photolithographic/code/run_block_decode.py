#!/usr/bin/env python3
"""Decode assigned pools, retaining global addresses and reference-blind search."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020';OUT=DATA/'results/block_decode'
REF=DATA/'sources/File3_synthesized.fasta'
COMMON={15:[(1,1.8),(1,2.55)],17:[(1,2.55),(1,3.5)],18:[(1,2.55),(1,3.5)],19:[(1,2.55),(1,3.5)],20:[(1,3.5),(1,5.0)]}

def job(mode,block,k,grid,reads,limit=200):
    folder=OUT/'decodes'/mode/f'block{block:02}';folder.mkdir(parents=True,exist_ok=True)
    prefix=folder/'file3';record=folder/f'k{k}.run.json'
    expected=[(c,float(r)) for c,r in grid]
    if record.exists():
        old=json.loads(record.read_text())
        if old['exit_code']==0 and old['grid']==[list(p) for p in expected]:
            print('SKIP',mode,block,k,flush=True);return
    env=os.environ.copy();env.update(ANTK_LAYOUT='internal',ANTK_BLOOM='1',ANTK_READ_LIMIT=str(limit),ANTK_GRID=','.join(f'{r:.2f}:{c}' for c,r in expected))
    cmd=[str(DATA/'bin/graph_block_rs'),str(reads),'196596',str(k),str(prefix),str(REF)]
    print('START',mode,'block',block,'k',k,'configurations',len(grid),flush=True);start=time.monotonic()
    with (folder/f'k{k}.log').open('w') as log:p=subprocess.run(['/usr/bin/time','-l']+cmd,env=env,stdout=log,stderr=log)
    result=dict(mode=mode,block=block,k=k,grid=expected,read_limit=limit,reads=str(reads),read_bytes=reads.stat().st_size,
                command=cmd,exit_code=p.returncode,seconds=time.monotonic()-start)
    if p.returncode==0:
        for c,r in expected:
            stem=folder/f'file3.k{k}.c{c}.r{r:.2f}'
            assert Path(str(stem)+'.strands').read_bytes()==Path(str(stem)+'.strict.strands').read_bytes()
            with Path(str(stem)+'.truth_candidates.tsv').open('rb') as handle:
                handle.seek(max(0,Path(str(stem)+'.truth_candidates.tsv').stat().st_size-256))
                assert handle.read().splitlines()[-1].startswith(b'196595\t')
    record.write_text(json.dumps(result,indent=2))
    print('FINISHED',mode,block,k,'exit',p.returncode,'seconds',round(result['seconds'],2),flush=True)
    if p.returncode:raise RuntimeError(result)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['controls','blocks','assigned-control'],required=True);ap.add_argument('--parallel',type=int,default=3)
    ap.add_argument('--include-unfiltered-raw',action='store_true',help='Optional expensive all-raw mixed-pool auxiliary controls')
    args=ap.parse_args()
    if args.phase=='controls':
        job('length60_control',0,19,[(1,3.5)],DATA/'processed/file3.fa.gz',60)
        prefix=OUT/'decodes/length60_control/block00/file3.k19.c1.r3.50'
        stats=json.loads(Path(str(prefix)+'.inner_stats.json').read_text());checks=[]
        for policy in stats['emitted']:
            a=Path(str(prefix)+f'.{policy}.strands');b=DATA/f'results/inner_rs/file3.k19.c1.r3.50.{policy}.strands'
            assert a.read_bytes()==b.read_bytes(),policy
            checks.append(dict(policy=policy,byte_identical=True,sha256=hashlib.sha256(a.read_bytes()).hexdigest()))
        (OUT/'read_limit_control.json').write_text(json.dumps(checks,indent=2))
        tasks=[('pooled_raw',0,k,grid,DATA/'raw/File3.fq.gz') for k,grid in COMMON.items()] if args.include_unfiltered_raw else []
    elif args.phase=='assigned-control':
        while not (OUT/'split_provenance.json').exists():time.sleep(20)
        pooled=OUT/'reads/assigned_pooled.fq.gz'
        with pooled.open('wb') as dest:
            for b in range(1,13):
                with (OUT/f'reads/broad/block{b:02}.fq.gz').open('rb') as src:shutil.copyfileobj(src,dest,1024*1024)
        stats=json.loads((OUT/'reads/split_counts.json').read_text())
        (OUT/'assigned_pooled.json').write_text(json.dumps(dict(reads=sum(r['broad_reads'] for r in stats['blocks']),
            path=str(pooled),bytes=pooled.stat().st_size,sha256=hashlib.sha256(pooled.read_bytes()).hexdigest(),
            construction='Binary concatenation of all 12 broad FASTQ gzip members; same sequences and qualities as the separate pools'),indent=2))
        tasks=[('assigned_pooled',0,k,grid,pooled) for k,grid in COMMON.items()]
    else:
        while not (OUT/'split_provenance.json').exists():time.sleep(20)
        plan=json.loads((DATA/'results/inner_rs/file3_candidate_path_bounds.json').read_text())
        tasks=[]
        for k in [19,18,17,20,15,13,21,25,31]:
            grid=sorted([(r['cutoff'],r['ratio']) for r in plan if r['k']==k],key=lambda p:(p[1],p[0]))
            tasks.append(('broad',12,k,grid,OUT/'reads/broad/block12.fq.gz'))
        for mode in ['strict','inclusive']:
            reads=OUT/'reads/strict/block12.fq.gz' if mode=='strict' else OUT/'reads/block12.inclusive.fq.gz'
            for k,grid in COMMON.items():tasks.append((mode,12,k,grid,reads))
        for b in [1,11,6,9,10,2,3,4,5,7,8]:
            for k,grid in COMMON.items():tasks.append(('broad',b,k,grid,OUT/f'reads/broad/block{b:02}.fq.gz'))
    with concurrent.futures.ThreadPoolExecutor(args.parallel) as pool:list(pool.map(lambda a:job(*a),tasks))
    print('DECODING PHASE COMPLETE',args.phase,flush=True)

if __name__=='__main__':main()
