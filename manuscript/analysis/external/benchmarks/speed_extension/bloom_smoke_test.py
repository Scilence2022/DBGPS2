"""Smoke-test both existing Bloom APIs through the new format caller."""
import json,os,subprocess
from pathlib import Path
R=Path(__file__).resolve().parent;S=R/'smoke'
configs=[('none',0,0),('B8',8,0),('B16',16,0),('F32',0,32),('F34',0,34),('B12_F34',12,34)];records=[]
for mode,k,c,ratio in [('fountain',21,0,6.),('nbt',29,0,1.9),('song',31,1,2.55)]:
 reads=S/('fountain.graph.fa' if mode=='fountain' else mode+'.fa');want=set((S/f'{mode}.graph.txt').read_text().splitlines())
 for name,qb,fb in configs:
  out=S/f'{mode}.{name}.bloom.txt';env=os.environ.copy();env.update(DBGPS_DETECT_ONLY='1',SPEED_QUERY_BITS=str(qb),SPEED_FILTER_SHIFT=str(fb));cmd=[str(R/'bin/format_graph_bloom'),mode,str(reads),str(R/'src/legacy-2base.map'),str(k),str(c),str(ratio),str(out),'0']
  subprocess.run(cmd,check=True,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);got=set(out.read_text().splitlines());assert got==want,(mode,name,len(got),len(want));records.append(dict(format=mode,config=name,exact=100,expected=100));print('BLOOM_TEST',mode,name,'100/100',flush=True)
(R/'bloom_smoke_validation.json').write_text(json.dumps(records,indent=2))
