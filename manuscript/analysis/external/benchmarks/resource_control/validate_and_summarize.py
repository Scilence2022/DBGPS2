"""Run only after all162 timed processes have finished."""
from pathlib import Path
import csv,hashlib,json,statistics,sys
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension';sys.path.insert(0,str(B));from run_bloom_matrix import CASES
DATASETS=['fountain','nbt_file16','song_ht176'];LABELS={'fountain':'DNA Fountain','nbt_file16':'Organick File 16','song_ht176':'HT176'}
def write_csv(path,rows):
 with path.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def main():
 rows=json.loads((R/'resource_results.json').read_text());assert len(rows)==162
 original={(z['dataset'],z['case']):z for z in json.loads((B/'matrix_results.json').read_text())};oldfile={z['output_sha256']:z for z in json.loads((B/'fountain_file_validation.json').read_text())};checks=[];summary=[];stages=[];allrecords=[]
 for workers,resource_mode in [(1,'serial'),(3,'three_workers')]:
  for dataset in DATASETS:
   baseline=statistics.mean(z['wall_seconds'] for z in rows if z['workers']==workers and z['dataset']==dataset and z['case']=='cluster_msa_vote')
   baseline_cpu=statistics.mean(z['cpu_seconds'] for z in rows if z['workers']==workers and z['dataset']==dataset and z['case']=='cluster_msa_vote')
   get=lambda case:[z for z in rows if z['workers']==workers and z['dataset']==dataset and z['case']==case]
   for case in CASES:
    group=get(case);assert len(group)==3 and {z['repeat'] for z in group}=={1,2,3};assert len({z['sha256'] for z in group})==1;expected=original[(dataset,case)];assert group[0]['sha256']==expected['sha256'];assert len({(z['emitted'],z['exact'],z['incorrect']) for z in group})==1
    for z in group:
     assert z['exit_code']==0 and z['cpu_seconds']>0 and z['wall_seconds']>0;data=Path(z['output_path']).read_bytes();assert hashlib.sha256(data).hexdigest()==z['sha256'] and len(data.splitlines())==z['emitted'];assert all(v>=0 for k,v in z['stages'].items() if k.endswith('_cpu_seconds'))
     allrecords.append({k:z.get(k) for k in ['resource_mode','workers','pipeline_threads','dataset','case','repeat','wall_seconds','cpu_seconds','user_seconds','system_seconds','native_wall_seconds','native_user_seconds','native_system_seconds','harness_user_seconds','harness_system_seconds','emitted','exact','incorrect','recovery','maximum_rss_bytes','reported_graph_rss_gb','k','cutoff','ratio','query_bits','filter_shift','query_bytes','sha256']})
     if case=='cluster_msa_vote':allrecords[-1]['pipeline_threads']=None
    z=group[0];mu=statistics.mean(x['wall_seconds'] for x in group);cpu=statistics.mean(x['cpu_seconds'] for x in group);rss=[x['reported_graph_rss_gb'] for x in group if x['reported_graph_rss_gb'] is not None]
    payload=None
    if dataset=='fountain':
     filecheck=oldfile[z['sha256']];assert filecheck['md5_match'] and filecheck['payload_mismatches_against_reencoded_archive']==0;payload=0
    summary.append(dict(workers=workers,pipeline_threads=z['pipeline_threads'] if case!='cluster_msa_vote' else None,resource_mode=resource_mode,dataset=dataset,label=LABELS[dataset],case=case,wall_mean=mu,wall_sd=statistics.stdev(x['wall_seconds'] for x in group),cpu_mean=cpu,cpu_sd=statistics.stdev(x['cpu_seconds'] for x in group),wall_speedup=baseline/mu,cpu_ratio_vs_msa=baseline_cpu/cpu,emitted=z['emitted'],exact=z['exact'],incorrect=z['incorrect'],recovery=z['recovery'],fountain_payload_inconsistent=payload,peak_graph_rss_gb=max(rss) if rss else None,k=z['k'],cutoff=z['cutoff'],ratio=z['ratio'],query_bits=z['query_bits'],filter_shift=z['filter_shift'],query_bytes=z.get('query_bytes',0),sha256=z['sha256']))
    checks.append(dict(workers=workers,dataset=dataset,case=case,repeats=3,identical_within_and_across_resource_settings=True,matches_original_eight_worker_output=True,sha256=z['sha256']))
    phase_names=['cluster','msa_vote','validation'] if case=='cluster_msa_vote' else ['count','bloom','seed','search']
    phase_wall=phase_cpu=0
    for phase in phase_names:
     wall_values=[x['stages'].get(phase+'_seconds',0) for x in group];cpu_values=[x['stages'].get(phase+'_cpu_seconds',0) for x in group];wm=statistics.mean(wall_values);cm=statistics.mean(cpu_values);phase_wall+=wm;phase_cpu+=cm
     stages.append(dict(workers=workers,resource_mode=resource_mode,dataset=dataset,case=case,stage=phase,wall_mean=wm,wall_sd=statistics.stdev(wall_values),cpu_mean=cm,cpu_sd=statistics.stdev(cpu_values)))
    stages.append(dict(workers=workers,resource_mode=resource_mode,dataset=dataset,case=case,stage='other',wall_mean=mu-phase_wall,wall_sd=None,cpu_mean=cpu-phase_cpu,cpu_sd=None))
   for case in ['dbgps2_B8','dbgps2_B16']:assert get(case)[0]['sha256']==get('dbgps2')[0]['sha256']
   for case in ['dbgps2_F34','dbgps2_B12_F34']:assert get(case)[0]['sha256']==get('dbgps2_F32')[0]['sha256']
 report=dict(completed_runs=162,all_exit_zero=True,all_repeats_identical=True,all_outputs_match_original_eight_worker_sets=True,query_filters_preserve_outputs=True,counting_filters_agree=True,fountain_file_validation='All new outputs are SHA256-identical to the original sets independently verified by complete archive MD5, payload re-encoding and RS re-encoding; no re-decoding is included in timing',checks=checks)
 (R/'resource_validation.json').write_text(json.dumps(report,indent=2));(R/'resource_summary.json').write_text(json.dumps(summary,indent=2));write_csv(R/'all_repeats_wall_and_cpu.csv',allrecords);write_csv(R/'all_configurations_summary.csv',summary);write_csv(R/'stage_wall_and_cpu.csv',stages)
 write_csv(R/'Table7_primary_methods.csv',[z for z in summary if z['case'] in ['cluster_msa_vote','dbgps_rule','dbgps2']]);write_csv(R/'Table8_primary_stage_timings.csv',[z for z in stages if z['case'] in ['cluster_msa_vote','dbgps_rule','dbgps2']]);print('VALIDATED_AND_SUMMARIZED',len(rows),len(summary),flush=True)
if __name__=='__main__':main()
