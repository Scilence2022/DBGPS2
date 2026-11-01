"""Check simulator accounting and profiler against actual native search."""
from pathlib import Path
import csv,json,subprocess
from simulate import *
R=Path(__file__).resolve().parent;D=R/'smoke';D.mkdir(exist_ok=True)
ref=designs(19,16);write_fasta(D/'ref.fa',ref.items());checks=[]
for length,depth,q,e,channel in [(200,10,1,0,'substitution'),(45,50,1,.01,'substitution'),(80,50,.05,.03,'substitution'),(25,30,1,0,'substitution'),(45,50,1,.03,'mixed')]:
 reads,meta=fragments(ref,19,length,depth,q);corrupted,errors=corrupt(reads,seed_of(19,length,depth,q,channel),e,channel);path=D/'reads.fa';write_fasta(path,corrupted)
 assert meta['retained_bases']==sum(len(s) for _,s in reads)
 assert errors['output_bases']==errors['source_bases']-errors['deletions']+errors['insertions']
 for k in [21,31]:
  prefix=D/f'L{length}_q{q}_e{e}_{channel}_k{k}'
  cmd=[str(R/'bin/boundary_profile'),'-x',str(D/'ref.fa'),'-a','1000','-b','1015','-k',str(k),'-t','1','-N','-v','1','-c','0,1,3','-r','1.25,2.55,4,0','-o',str(prefix),str(path)]
  p=subprocess.run(cmd,capture_output=True,text=True);assert p.returncode==0,(cmd,p.stderr)
  rows=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(rows)==12 and all(None not in z for z in rows)
  assert all(int(z['native_checks'])==16 for z in rows)
  if length==200 and e==0:assert all(int(z['correct'])==16 for z in rows)
  if length==25 and k==31:assert all(int(z['ceil_cov'])==0 and int(z['correct'])==0 for z in rows)
  checks.append(dict(length=length,k=k,channel=channel,native_search_checks=192))
# Integer truncation differs from the conservative continuous ratio condition.
a,b,r=6,2,2.55;assert b>=int(a/r) and b<=int(a*r) and a/b>r
(R/'smoke_validation.json').write_text(json.dumps(dict(passed=True,independent_native_search_comparisons=sum(x['native_search_checks'] for x in checks),checks=checks),indent=2));print('SMOKE_OK',len(checks),sum(x['native_search_checks'] for x in checks))
