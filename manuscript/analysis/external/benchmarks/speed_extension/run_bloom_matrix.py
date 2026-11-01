"""Expanded serial matrix: conventional pipeline versus graph configurations.
Reuses complete measurements of unchanged no-Bloom comparators; no partial
run contributes a timing. Every newly added configuration is timed three times.
"""
import hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
import run_benchmark as base
R=Path(__file__).resolve().parent
CASES=['cluster_msa_vote','dbgps_rule','dbgps2','dbgps2_alt_c','dbgps2_B8','dbgps2_B16','dbgps2_F32','dbgps2_F34','dbgps2_B12_F34']
SETTINGS={'dbgps2_B8':(8,0),'dbgps2_B16':(16,0),'dbgps2_F32':(0,32),'dbgps2_F34':(0,34),'dbgps2_B12_F34':(12,34)}
def run(label,case,rep,chosen):
 cfg=base.DATA[label]
 if case in ['cluster_msa_vote','dbgps_rule','dbgps2']:
  cut=chosen[label]['dbgps_cutoff'] if case=='dbgps_rule' else cfg['c'];row=base.execute(label,case,rep,cut);row=dict(row,query_bits=0,filter_shift=0,case=case);return row
 qb,fb=SETTINGS.get(case,(0,0));cut=(1 if cfg['c']==0 else chosen[label]['dbgps_cutoff']) if case=='dbgps2_alt_c' else cfg['c']
 prefix=R/'runs'/f'{label}_{case}_c{cut}_rep{rep}';result=Path(str(prefix)+'.result.json')
 if result.exists():return json.loads(result.read_text())
 env=base.ENV.copy();env.update(SPEED_QUERY_BITS=str(qb),SPEED_FILTER_SHIFT=str(fb));raw=Path(str(prefix)+'.raw.txt')
 cmd=[str(R/'bin/format_graph_bloom'),cfg['mode'],str(cfg.get('graph',cfg['reads'])),str(R/'src/legacy-2base.map'),str(cfg['k']),str(cut),str(cfg['r']),str(raw),'0']
 Path(str(prefix)+'.load_before.txt').write_text(subprocess.check_output(['ps','-Ao','pid,pcpu,comm','-r'],text=True));print('START',prefix.name,flush=True);start=time.perf_counter()
 with open(str(prefix)+'.stdout','w') as out,open(str(prefix)+'.stderr','w') as err:cp=subprocess.run(['/usr/bin/time','-l']+cmd,stdout=out,stderr=err,env=env)
 assert cp.returncode==0,(prefix,cp.returncode)
 dst=Path(str(prefix)+'.strands.txt');conflicts=base.canonicalize(raw,dst,cfg['mode']);elapsed=time.perf_counter()-start
 log=Path(str(prefix)+'.stderr').read_text();stages={x.lower()+'_seconds':float(re.search(x+r'.*?wall=([\d.]+)',log)[1]) for x in ['COUNT','BLOOM','SEED','SEARCH']}
 row=dict(dataset=label,method=case,case=case,repeat=rep,pilot=False,screen_modulus=1,wall_seconds=elapsed,cutoff=cut,ratio=cfg['r'],k=cfg['k'],query_bits=qb,filter_shift=fb,query_bytes=int(re.search(r'BLOOM bits=\d+ bytes=(\d+)',log)[1]),command=cmd,environment_options={k:env[k] for k in ['DBGPS_DETECT_ONLY','SPEED_QUERY_BITS','SPEED_FILTER_SHIFT']},exit_code=cp.returncode,conflicts=conflicts,stages=stages,**base.score(dst,cfg));result.write_text(json.dumps(row,indent=2));print('DONE',prefix.name,round(elapsed,2),'seconds',row['emitted'],row['exact'],row['incorrect'],flush=True);return row

def enrich(row):
 row=dict(row)
 method=row['method'];c=row['cutoff'] if row['cutoff'] is not None else base.DATA[row['dataset']]['c'];prefix=R/'runs'/f'{row["dataset"]}_{method}_c{c}_rep{row["repeat"]}'
 log=Path(str(prefix)+'.stderr').read_text();match=re.search(r'(\d+)\s+maximum resident set size',log)
 row['maximum_rss_bytes']=int(match[1]) if match else None
 row['reported_graph_rss_gb']=row['maximum_rss_bytes']/1e9 if match and method!='cluster_msa_vote' else None
 row['rss_scope']='single graph process and its threads' if method!='cluster_msa_vote' else 'time(1) process maximum; not summed over concurrent workers, excluded from memory comparison'
 row['output_path']=str(prefix)+'.strands.txt';row['result_path']=str(prefix)+'.result.json'
 return row

def main():
 chosen=json.loads((R/'selected_settings.json').read_text());records=[]
 metadata=dict(cases=CASES,repeats=3,datasets=list(base.DATA),host='same Apple M3 Max workstation',method_order='cyclic order across cases and datasets; completed unchanged comparator measurements retained from the documented earlier rotated-order block',serial=True,unmodified_completed_measurements_reused=True,comparison='elapsed full reconstruction after shared preprocessing; quality and resource use reported at every point',filters={'B':[8,16],'F':[32,34],'combined':{'B':12,'F':34}},parameter_control='same k, ratio and primary cutoff across all Bloom modes within a dataset; extra no-Bloom point uses c1 for Fountain/Organick, and the selected DBGPS cutoff for HT176, retaining the DBGPS2 ratio',rss='native graph process RSS only; MSA process-tree memory not inferred from time(1)',core='same saved counter, Bloom implementations, bounded traversal and format validators; a separate caller enables existing Bloom APIs; production sources unchanged')
 (R/'matrix_metadata.json').write_text(json.dumps(metadata,indent=2))
 for rep in range(1,4):
  for di,label in enumerate(base.DATA):
   shift=(rep-1+di*3)%len(CASES);order=CASES[shift:]+CASES[:shift]
   for case in order:
    row=enrich(run(label,case,rep,chosen));records.append(row);(R/'matrix_results.json').write_text(json.dumps(records,indent=2))
 print('MATRIX_ALL_DONE',len(records),flush=True)
if __name__=='__main__':main()
