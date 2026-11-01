#!/usr/bin/env python3
"""Run the user's DBGPS-analyzer on all File 3 blocks and read variants."""
import argparse
import concurrent.futures
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/analyzer_block12'
TOOLS=DATA/'bin/analyzer_tools'
REPO=Path('/Users/song/.codex/worktrees/dabb/DBGPS-analyzer')
REF=DATA/'sources/File3_synthesized.fasta'
KS=[13,15,17,18,19,20,21,25,31]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--parallel',type=int,default=3);args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    lines=REF.read_text().splitlines();assert len(lines)==196596*2
    for b in range(12):(OUT/f'block{b+1:02}.fasta').write_text('\n'.join(lines[b*32766:(b+1)*32766])+'\n')
    plan=[dict(dataset='processed',k=19,limit=60,bloom=False),
          dict(dataset='raw',k=19,limit=200,bloom=True),
          dict(dataset='processed',k=20,limit=60,bloom=True),
          dict(dataset='processed',k=17,limit=60,bloom=True),
          dict(dataset='processed',k=18,limit=60,bloom=True),
          dict(dataset='raw',k=19,limit=60,bloom=True)]
    plan += [dict(dataset='processed',k=k,limit=60,bloom=True) for k in KS if k not in [17,18,19,20]]
    plan += [dict(dataset='raw',k=k,limit=200,bloom=True) for k in [17,20]]
    (OUT/'plan.json').write_text(json.dumps(plan,indent=2))
    def run(row):
        name=f"{row['dataset']}.k{row['k']}.L{row['limit']}";prefix=OUT/name
        record=OUT/(name+'.run.json')
        if record.exists() and json.loads(record.read_text()).get('exit_code')==0:
            print('SKIP',name,flush=True);return
        reads=DATA/('processed/file3.fa.gz' if row['dataset']=='processed' else 'raw/File3.fq.gz')
        cmd=[str(TOOLS/'DBGPS-analyzer'),'-k',str(row['k']),'-t','8','-L',str(row['limit']),'-c','0','-C','3','-r','0','-o',str(prefix)]
        if row['bloom']:cmd+=['--bloom']
        cmd += [str(REF),str(reads)]
        print('START',name,'reference-only' if row['bloom'] else 'full read index',flush=True)
        start=time.monotonic()
        with (OUT/(name+'.stdout.log')).open('w') as stdout,(OUT/(name+'.stderr.log')).open('w') as stderr:
            p=subprocess.run(['/usr/bin/time','-l']+cmd,stdout=stdout,stderr=stderr)
        result=dict(**row,command=cmd,exit_code=p.returncode,seconds=time.monotonic()-start,
                    reads=str(reads),read_bytes=reads.stat().st_size)
        record.write_text(json.dumps(result,indent=2))
        print('FINISHED',name,'exit',p.returncode,'seconds',round(result['seconds'],2),flush=True)
        if p.returncode:raise RuntimeError(name)
        with (OUT/(name+'.cov_details')).open() as handle:assert sum(1 for _ in handle)==196596
    with concurrent.futures.ThreadPoolExecutor(args.parallel) as pool:list(pool.map(run,plan))
    def links(k):
        results=[]
        for b in range(13):
            ref=REF if b==0 else OUT/f'block{b:02}.fasta'
            p=subprocess.run([str(TOOLS/'DBGPS-links'),'-k',str(k),'-p','0',str(ref)],capture_output=True,text=True,check=True)
            m=re.search(r'Total cross links (\d+)',p.stdout);assert m
            results.append(dict(block=b,cross_links=int(m[1]),strands=196596 if not b else 16383))
        (OUT/f'links.k{k}.json').write_text(json.dumps(results,indent=2));print('LINKS',k,'complete',flush=True)
    with concurrent.futures.ThreadPoolExecutor(3) as pool:list(pool.map(links,KS))
    sources=['DBGPS-analyzer.c','DBGPS-links.c','dbgps_core.h','khashl.h','kseq.h','kthread.c','vendor/bfc/bbf.c','vendor/bfc/bbf.h']
    provenance=dict(repository=str(REPO),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        source_sha256={s:hashlib.sha256((REPO/s).read_bytes()).hexdigest() for s in sources},
        binary_sha256={s:hashlib.sha256((TOOLS/s).read_bytes()).hexdigest() for s in ['DBGPS-analyzer','DBGPS-links']},
        reference=str(REF),reference_sha256=hashlib.sha256(REF.read_bytes()).hexdigest(),
        build_flags='cc -g -Wall -Wextra -O3; -lz -lpthread',primer_trim_front=0,primer_trim_back=0,
        count_cap=16383,reference_only_contract='Exact target counts; global off-target noise is unavailable for Bloom runs')
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2))
    print('ALL ANALYZER RUNS COMPLETE',flush=True)

if __name__=='__main__':main()
