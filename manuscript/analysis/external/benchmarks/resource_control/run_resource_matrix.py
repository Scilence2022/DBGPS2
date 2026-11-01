"""New complete serial/three-worker timings; previous 8-worker outputs are quality controls only."""
from pathlib import Path
import argparse,hashlib,json,os,platform,re,resource,subprocess,sys,time
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension';sys.path.insert(0,str(B))
import run_benchmark as base
from run_bloom_matrix import CASES,SETTINGS

def parent_cpu():
 u=resource.getrusage(resource.RUSAGE_SELF);return u.ru_utime,u.ru_stime

def atomic_json(path,value):
 tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2));tmp.replace(path)

def execute(dataset,case,repeat,workers,chosen,expected):
 cfg=base.DATA[dataset];mode='serial' if workers==1 else 'three_workers';folder=R/mode/'runs';folder.mkdir(parents=True,exist_ok=True)
 prefix=folder/f'{dataset}_{case}_rep{repeat}';result=Path(str(prefix)+'.result.json')
 if result.exists():
  row=json.loads(result.read_text());assert row['exit_code']==0 and row['workers']==workers and row['sha256']==expected['sha256'];return row
 env=base.ENV.copy();env.update(BENCHMARK_COUNT_WORKERS=str(workers),BENCHMARK_PIPELINE_THREADS='1' if workers==1 else '3')
 qb,fb=SETTINGS.get(case,(0,0));cut=chosen[dataset]['dbgps_cutoff'] if case=='dbgps_rule' else (1 if cfg['c']==0 else chosen[dataset]['dbgps_cutoff']) if case=='dbgps2_alt_c' else cfg['c'];ratio=1000000. if case=='dbgps_rule' else cfg['r']
 env.update(SPEED_QUERY_BITS=str(qb),SPEED_FILTER_SHIFT=str(fb))
 if case=='cluster_msa_vote':
  command=[sys.executable,str(R/'msa_pipeline.py'),cfg['mode'],str(cfg['reads']),str(prefix),'--workers',str(workers)];raw=Path(str(prefix)+'.valid.txt')
 else:
  raw=Path(str(prefix)+'.raw.txt');binary='format_graph' if case in ['dbgps_rule','dbgps2'] else 'format_graph_bloom';command=[str(R/'bin'/binary),cfg['mode'],str(cfg.get('graph',cfg['reads'])),str(R/'src/legacy-2base.map'),str(cfg['k']),str(cut),str(ratio),str(raw),'0']
 Path(str(prefix)+'.load_before.txt').write_text(subprocess.check_output(['ps','-Ao','pid,pcpu,comm','-r'],text=True))
 print('START',mode,prefix.name,time.strftime('%Y-%m-%dT%H:%M:%S%z'),flush=True);before=parent_cpu();start=time.perf_counter()
 with open(str(prefix)+'.stdout','w') as out,open(str(prefix)+'.stderr','w') as err:cp=subprocess.run(['/usr/bin/time','-l']+command,stdout=out,stderr=err,env=env)
 assert cp.returncode==0,(prefix,cp.returncode)
 dst=Path(str(prefix)+'.strands.txt');conflicts=base.canonicalize(raw,dst,cfg['mode']);elapsed=time.perf_counter()-start;after=parent_cpu()
 log=Path(str(prefix)+'.stderr').read_text();time_match=re.search(r'([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys',log);assert time_match,prefix;native_wall,user,system=map(float,time_match.groups());parent_user=after[0]-before[0];parent_system=after[1]-before[1]
 if case=='cluster_msa_vote':stages=json.loads(Path(str(prefix)+'.stages.json').read_text())
 else:
  stages={}
  for tag in ['COUNT','BLOOM','SEED','SEARCH']:
   match=re.search(r'^'+tag+r'.*?wall=([\d.]+) cpu=([\d.]+)',log,re.M)
   if match:stages[tag.lower()+'_seconds']=float(match[1]);stages[tag.lower()+'_cpu_seconds']=float(match[2])
  stages['reconstruction_process_seconds']=native_wall
 rss_match=re.search(r'(\d+)\s+maximum resident set size',log);query_match=re.search(r'BLOOM bits=\d+ bytes=(\d+)',log)
 score=base.score(dst,cfg);assert score['sha256']==expected['sha256'],(mode,dataset,case,'output differs from previous verified configuration')
 row=dict(resource_mode=mode,workers=workers,pipeline_threads=1 if workers==1 else 3,case=case,method=case,dataset=dataset,repeat=repeat,wall_seconds=elapsed,native_wall_seconds=native_wall,native_user_seconds=user,native_system_seconds=system,harness_user_seconds=parent_user,harness_system_seconds=parent_system,user_seconds=user+parent_user,system_seconds=system+parent_system,cpu_seconds=user+system+parent_user+parent_system,exit_code=0,command=command,environment_options={k:env[k] for k in ['DBGPS_DETECT_ONLY','SPEED_QUERY_BITS','SPEED_FILTER_SHIFT','BENCHMARK_COUNT_WORKERS','BENCHMARK_PIPELINE_THREADS']},cutoff=cut if case!='cluster_msa_vote' else None,ratio=ratio if case!='cluster_msa_vote' else None,k=cfg['k'] if case!='cluster_msa_vote' else None,query_bits=qb,filter_shift=fb,query_bytes=int(query_match[1]) if query_match else 0,maximum_rss_bytes=int(rss_match[1]) if rss_match else None,reported_graph_rss_gb=int(rss_match[1])/1e9 if rss_match and case!='cluster_msa_vote' else None,stages=stages,conflicts=conflicts,output_path=str(dst),result_path=str(result),matches_previous_eight_worker_output=True,**score)
 atomic_json(result,row);print('DONE',mode,prefix.name,round(elapsed,2),'wall_s',round(row['cpu_seconds'],2),'cpu_s',score['emitted'],flush=True);return row

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--only-workers',type=int,choices=[1,3]);a=ap.parse_args();chosen=json.loads((B/'selected_settings.json').read_text());old=json.loads((B/'matrix_results.json').read_text());expected={(x['dataset'],x['case']):x for x in old};rows=[];workersets=[1,3] if a.only_workers is None else [a.only_workers]
 metadata=dict(datasets=list(base.DATA),cases=CASES,repeats=3,worker_settings=workersets,target_runs=3*3*9*len(workersets),inputs='same complete retained molecules as previous 81-run matrix; all parameters and validators unchanged',serial_protocol='count_workers=1,pipeline_threads=1; no overlapping graph pipeline stages; MSA direct sequential map with no Python process pool',three_worker_protocol='count_workers=3,pipeline_threads=3; MSA workers=3; primary worker counts matched, not a strict three-core limit',timing='start through native reconstruction, format validation, reference-free conflict rejection and sorted output writing; reference truth scoring and outer decoding follow the timer',cpu='time(1) user+system plus calling-harness self CPU over the same interval; getrusage stage times include graph threads or completed MSA worker/child usage',threading='numerical library parallelism capped at one; MUSCLE3 calls single-threaded',order='all measured jobs serial; cyclic method order and alternating resource-setting order within each dataset/repeat; no previous timings reused',host=platform.platform(),dedicated_host=False,cpu_affinity=False,background_applications='present; pre-run snapshots retained',reference_evaluation='same design scoring as original; original outputs used only after timing to check scheduler invariance')
 atomic_json(R/'resource_metadata.json',metadata)
 for rep in range(1,4):
  for di,dataset in enumerate(base.DATA):
   shift=(rep-1+di*3)%len(CASES);cases=CASES[shift:]+CASES[:shift]
   for ci,case in enumerate(cases):
    order=workersets if (rep+di+ci)%2 else workersets[::-1]
    for workers in order:
     row=execute(dataset,case,rep,workers,chosen,expected[(dataset,case)]);rows.append(row);atomic_json(R/'resource_results.json',rows)
 print('RESOURCE_MATRIX_ALL_DONE',len(rows),flush=True)
if __name__=='__main__':main()
