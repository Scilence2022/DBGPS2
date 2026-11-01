#!/usr/bin/env python3
"""Reference-assisted partitioning of all raw File 3 reads; stream SAM only."""
import hashlib
import json
import subprocess
import time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'input/antkowiak2020'
OUT=DATA/'results/block_decode'

def main():
    target=OUT/'reads'
    for sub in ['broad','strict']:(target/sub).mkdir(parents=True,exist_ok=True)
    reference=DATA/'sources/File3_synthesized.fasta';reads=DATA/'raw/File3.fq.gz'
    mapper=['minimap2','-a','-x','sr','-k','11','-w','2','-n','2','-m','15','-s','20',
            '-A','2','-B','3','-O','4,12','-E','2,1','-N','5','-p','0.8','--secondary=yes','--eqx','-t','12',str(reference),str(reads)]
    router=[str(DATA/'bin/split_aligned_blocks'),str(target)]
    start=time.monotonic()
    with (OUT/'mapping.log').open('w') as mlog,(OUT/'splitting.log').open('w') as rlog:
        m=subprocess.Popen(['/usr/bin/time','-l']+mapper,stdout=subprocess.PIPE,stderr=mlog)
        r=subprocess.Popen(['/usr/bin/time','-l']+router,stdin=m.stdout,stderr=rlog)
        m.stdout.close()
        router_exit=r.wait()
        if router_exit:m.terminate()
        mapper_exit=m.wait()
    assert mapper_exit==router_exit==0,(mapper_exit,router_exit)
    stats=json.loads((target/'split_counts.json').read_text())
    assert stats['total_reads']==195619515
    assert sum(b['broad_reads'] for b in stats['blocks'])+stats['ambiguous_reads']+stats['unassigned_reads']==stats['total_reads']
    assert all(b['strict_reads']<=b['broad_reads'] for b in stats['blocks'])
    provenance=dict(mapper_command=mapper,splitter_command=router,mapper_version=subprocess.check_output(['minimap2','--version'],text=True).strip(),
        seconds=time.monotonic()-start,mapper_exit=mapper_exit,router_exit=router_exit,
        reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),
        router_source_sha256=hashlib.sha256((ROOT/'benchmarks/antkowiak2020/split_aligned_blocks.cpp').read_bytes()).hexdigest(),
        router_binary_sha256=hashlib.sha256((DATA/'bin/split_aligned_blocks').read_bytes()).hexdigest(),
        reference_assistance='Read assignment uses the known full design payload; this is not a reference-free end-to-end reconstruction experiment',
        broad='Reported alignments: matches>=20, query-aligned bases>=30, payload overlap>=20, AS>=25, match fraction>=0.6; every qualifying hit with AS>=0.8 of best qualifying score must name the same block',
        strict='Subset of broad: primary MAPQ>=10, matches>=24, query-aligned bases>=35, payload overlap>=25, AS>=40, match fraction>=0.7',
        ambiguous='Qualifying near-best reported alignments name multiple blocks; at most five secondary alignments are reported',
        inclusive='Block 12 broad reads plus ambiguous reads with a qualifying near-best block 12 alignment',
        sequence_handling='Original full sequences and qualities, reverse-complement/reversal undone on reverse-strand SAM primaries',
        outputs={str(p.relative_to(OUT)):p.stat().st_size for p in target.rglob('*.fq.gz')},counts=stats)
    (OUT/'split_provenance.json').write_text(json.dumps(provenance,indent=2))
    print(json.dumps(stats,indent=2));print('FULL RAW READ PARTITION COMPLETE',flush=True)

if __name__=='__main__':main()
