#!/usr/bin/env python3
"""Offline compilation and perfect-input controls for the packaged C adapters."""
from pathlib import Path
import csv, json, subprocess, tempfile

ROOT=Path(__file__).resolve().parent
EXP=ROOT/'experiments/nbt4079'
COMP=str.maketrans('ACGT','TGCA')
def rc(s):return s.translate(COMP)[::-1]
seqs=[]
with (ROOT/'fixtures/perfect.fa').open() as h:
    for line in h:
        if not line.startswith('>') and line.strip() not in seqs:seqs.append(line.strip())
assert len(seqs)==20
truth={s[20:130] for s in seqs}
results={}
with tempfile.TemporaryDirectory(prefix='dbgps2-package-') as tmp:
    tmp=Path(tmp)
    seedfile=tmp/'both.seeds.tsv'
    seedfile.write_text(''.join(f'{direction}\t{oriented[3:34]}\t3\n'
        for s in seqs for direction,oriented in [('F',s),('R',rc(s))]))
    refs=tmp/'refs.fa';refs.write_text(''.join(f'>ref_{i}\n{s}\n' for i,s in enumerate(seqs)))
    gap=tmp/'gap.fa'
    altered=[s[:75]+('C' if s[75]!='C' else 'A')+s[76:] for s in seqs]
    gap.write_text(''.join(f'>gap_{i}_{j}\n{s}\n' for i,s in enumerate(altered) for j in range(3)))
    source=EXP/'optimization/scripts'
    for mode in ['native','inner','bound']:
        binary=tmp/mode
        cmd=['cc','-O3','-Wall','-I.']
        if mode=='inner':cmd+=['-DUSE_INNER']
        cmd+=['-o',str(binary),str(source/('graph_bound.c' if mode=='bound' else 'search_core.c')),
              'dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c','-lz','-lpthread']
        subprocess.run(cmd,cwd=ROOT,check=True)
        if mode!='bound':
            for k in [21,23,25,27,29,31]:
                output=tmp/f'{mode}.{k}.tsv'
                process=subprocess.run([str(binary),str(ROOT/'fixtures/perfect.fa'),str(k),'1'],
                    input=f'control 1 2.55 0 {seedfile} {output}\nQUIT\n',text=True,capture_output=True,check=True)
                with output.open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
                assert len(rows)==40 and all(r['payload'] in truth and len(r['payload'])==110 for r in rows)
                assert {r['payload'] for r in rows if r['direction']=='F'}==truth
                assert {r['payload'] for r in rows if r['direction']=='R'}==truth
                results[f'{mode}_k{k}']={'queries':40,'exact':40}
        else:
            for k in [29,31]:
                for name,reads,expected in [('perfect',ROOT/'fixtures/perfect.fa',20),('systematic_gap',gap,0)]:
                    process=subprocess.run([str(binary),str(reads),str(seedfile),str(k),str(refs)],
                                           text=True,capture_output=True,check=True)
                    result=json.loads(process.stdout)
                    assert result['bounds'][0]['either_direction']==expected
                    results[f'bound_{name}_k{k}']={'expected':expected,'observed':result['bounds'][0]['either_direction']}
(ROOT/'ADAPTER_VALIDATION.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
