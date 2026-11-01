#!/usr/bin/env python3
"""Portable staged reproduction. Run with Python >=3.10 from the extracted folder."""
import argparse, gzip, hashlib, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if any(c.isspace() for c in str(ROOT)):
    raise SystemExit('Extract into a path without whitespace; the historical C driver protocol uses whitespace-delimited paths.')
EXP=ROOT/'experiments/nbt4079'
for rel in ['data','sources','tmp','results','optimization/tmp','optimization/results']:
    (EXP/rel).mkdir(parents=True,exist_ok=True)
def run(rel,*args):
    subprocess.run([sys.executable,str(EXP/rel),*args],cwd=ROOT,check=True)
def references():
    n=0
    with gzip.open(EXP/'data/id20.refs.txt.gz','rt') as source,(EXP/'data/references.fa').open('w') as dest:
        for n,line in enumerate(source,1):
            assert line.startswith("5'-") and line.rstrip().endswith("-3'")
            seq=line.strip()[3:-3]
            assert len(seq)==150 and set(seq)<=set('ACGT')
            dest.write(f'>ref_{n-1}\n{seq}\n')
    assert n==607150
def smoke():
    binary=EXP/'tmp/smoke_adapter'
    cmd=['cc','-O3','-Wall','-I.','-o',str(binary),str(EXP/'scripts/dbgps2_adapter.c'),
         'dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c','-lz','-lpthread']
    subprocess.run(cmd,cwd=ROOT,check=True)
    bodies={}
    with (ROOT/'fixtures/perfect.fa').open() as h:
        for line in h:
            if not line.startswith('>'):
                body=line.strip()[20:130];bodies[body[:14]]=body
    result={}
    for stem,kind in [('perfect','perfect'),('noisy_control','noisy')]:
        prefix=EXP/'tmp'/('smoke_'+kind)
        subprocess.run([str(binary),str(ROOT/'fixtures'/f'{stem}.fa'),str(ROOT/'fixtures'/f'{stem}.seeds.tsv'),str(prefix),'1'],cwd=ROOT,check=True)
        import csv
        configs=['c0_r2.55','c1_r2.55'] if kind=='perfect' else ['c1_r2.55']
        for config in configs:
            with Path(str(prefix)+'.'+config+'.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
            exact=sum(r['payload']==bodies.get(r['seed'][-14:]) for r in rows)
            assert exact==len(rows)==20
            result[kind+'_'+config]={'exact':exact,'returned':len(rows)}
    (ROOT/'SMOKE_RESULT.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
def check():
    snap=json.loads((ROOT/'SOURCE_SNAPSHOT.json').read_text())
    for name,sha in snap['source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    print('Core source snapshot checksums match.')
stages=['check','download','count','preprocess','baseline','seeds','validate','optimize','inner','refine','aggregate','bounds','figures']
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--stage',choices=stages+['all','smoke'],default='check')
args=parser.parse_args()
for stage in stages if args.stage=='all' else [args.stage]:
    print('STAGE',stage,flush=True)
    if stage=='check':check()
    elif stage=='smoke':check();smoke()
    elif stage=='download':run('scripts/download_data.py')
    elif stage=='count':run('scripts/count_stream.py');references()
    elif stage=='preprocess':run('scripts/prepare_reads.py','--total','15719903')
    elif stage=='baseline':run('scripts/run_experiment.py');run('scripts/export_fasta.py')
    elif stage=='seeds':run('optimization/scripts/prepare_seeds.py')
    elif stage=='validate':run('optimization/scripts/validate_inner.py')
    elif stage=='optimize':run('optimization/scripts/optimize.py')
    elif stage=='inner':run('optimization/scripts/optimize_inner.py')
    elif stage=='refine':run('optimization/scripts/refine.py')
    elif stage=='aggregate':run('optimization/scripts/aggregate.py','--key-length','13','--include-refinement')
    elif stage=='bounds':run('optimization/scripts/compute_bounds.py')
    elif stage=='figures':run('scripts/plot_results.py');run('optimization/scripts/report.py')
