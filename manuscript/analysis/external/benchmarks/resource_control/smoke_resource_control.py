from pathlib import Path
import os,subprocess,sys,json
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension';S=B/'smoke';O=R/'smoke';O.mkdir(exist_ok=True);records=[]
configs=[('none',0,0),('B8',8,0),('B16',16,0),('F32',0,32),('F34',0,34),('B12_F34',12,34)]
for workers in [1,3]:
 env=os.environ.copy();env.update(DBGPS_DETECT_ONLY='1',BENCHMARK_COUNT_WORKERS=str(workers),BENCHMARK_PIPELINE_THREADS='1' if workers==1 else '3',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 for mode,k,c,ratio in [('fountain',21,0,6.),('nbt',29,0,1.9),('song',31,1,2.55)]:
  reads=S/('fountain.graph.fa' if mode=='fountain' else mode+'.fa');want=set((S/f'{mode}.graph.txt').read_text().splitlines())
  for name,qb,fb in configs:
   out=O/f'{mode}.w{workers}.{name}.txt';env.update(SPEED_QUERY_BITS=str(qb),SPEED_FILTER_SHIFT=str(fb));cmd=[str(R/'bin/format_graph_bloom'),mode,str(reads),str(R/'src/legacy-2base.map'),str(k),str(c),str(ratio),str(out),'0'];p=subprocess.run(cmd,check=True,env=env,capture_output=True,text=True);got=set(out.read_text().splitlines());assert got==want,(mode,name,workers,len(got),len(want));assert f'count_workers={workers} pipeline_threads={1 if workers==1 else 3}' in p.stderr;records.append(dict(format=mode,workers=workers,method=name,exact=len(got)));print('GRAPH_SMOKE',mode,workers,name,len(got),flush=True)
  prefix=O/f'{mode}.w{workers}.msa';cmd=[sys.executable,str(R/'msa_pipeline.py'),mode,str(S/(mode+'.fa')),str(prefix),'--workers',str(workers)];p=subprocess.run(cmd,check=True,env=env,capture_output=True,text=True);got=set(Path(str(prefix)+'.valid.txt').read_text().splitlines());assert got==want,(mode,workers,len(got),len(want));record=json.loads(Path(str(prefix)+'.stages.json').read_text());assert record['pipeline_cpu_seconds']>0;records.append(dict(format=mode,workers=workers,method='msa',exact=len(got),cpu_seconds=record['pipeline_cpu_seconds']));print('MSA_SMOKE',mode,workers,len(got),flush=True)
(R/'smoke_validation.json').write_text(json.dumps(records,indent=2));print('RESOURCE_SMOKES_COMPLETE',len(records),flush=True)
