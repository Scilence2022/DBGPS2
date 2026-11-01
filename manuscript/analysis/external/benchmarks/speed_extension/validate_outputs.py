"""Independent completed-run audit and untimed DNA Fountain file check."""
from pathlib import Path
import collections,hashlib,json,os,subprocess,sys
R=Path(__file__).resolve().parent
rows=json.loads((R/'results.json').read_text());assert len(rows)==27
checks=[]
for (dataset,method),group in __import__('itertools').groupby(sorted(rows,key=lambda x:(x['dataset'],x['method'])),key=lambda x:(x['dataset'],x['method'])):
 group=list(group);assert len(group)==3 and {x['repeat'] for x in group}=={1,2,3};assert all(x['exit_code']==0 and x['wall_seconds']>0 for x in group);assert len({x['sha256'] for x in group})==1
 for row in group:
  prefix=R/'runs'/f'{dataset}_{method}_c{row["cutoff"] if row["cutoff"] is not None else 0}_rep{row["repeat"]}'
  # conventional c suffix uses DBGPS2 c (one for HT176), although cutoff is not a method parameter
  if method=='cluster_msa_vote' and dataset=='song_ht176':prefix=R/'runs'/f'{dataset}_{method}_c1_rep{row["repeat"]}'
  p=Path(str(prefix)+'.strands.txt');b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==row['sha256'];assert len(b.splitlines())==row['emitted']
 checks.append(dict(dataset=dataset,method=method,repeats=3,identical=True,sha256=group[0]['sha256']))
# Expected hash is from the DNA Fountain published archive, already checked in full-raw-read analysis.
expected='8651e90d3a013178b816b63fdbb94b9b';filechecks=[]
for method in ['cluster_msa_vote','dbgps_rule','dbgps2']:
 row=next(x for x in rows if x['dataset']=='fountain' and x['method']==method and x['repeat']==1);c=row['cutoff'] if row['cutoff'] is not None else 0;p=R/'runs'/f'fountain_{method}_c{c}_rep1.strands.txt';out=R/'runs'/f'fountain_{method}.decoded.bin'
 cmd=[sys.executable,str(R/'src/decode_fountain.py'),str(p),'--no_screen','--out',str(out)]
 env=os.environ.copy();env['PYTHONPATH']='/tmp/dbgps2_speed_extension_dependencies'
 cp=subprocess.run(cmd,capture_output=True,text=True,env=env)
 (R/'runs'/f'fountain_{method}.decode.log').write_text(cp.stdout+cp.stderr)
 md5=hashlib.md5(out.read_bytes()).hexdigest() if out.exists() else None
 filechecks.append(dict(method=method,command=cmd,exit_code=cp.returncode,bytes=out.stat().st_size if out.exists() else None,md5=md5,expected_md5=expected,md5_match=md5==expected,decoding_in_timing=False))
 print(method,md5,md5==expected,flush=True)
# Re-encode each observed seed against the checksum-verified archive. This
# validates payload consistency without asserting original seed-list membership.
if all(x['md5_match'] for x in filechecks):
 sys.path.insert(0,str(R/'src'))
 from fountain_common import PRNG,dna_to_bytes
 blob=(R/'runs/fountain_dbgps2.decoded.bin').read_bytes();chunks=[blob[i:i+32] for i in range(0,len(blob),32)];prng=PRNG(67088,.001,.025);expected_by_seed={}
 for check in filechecks:
  row=next(x for x in rows if x['dataset']=='fountain' and x['method']==check['method'] and x['repeat']==1);c=row['cutoff'] if row['cutoff'] is not None else 0;path=R/'runs'/f'fountain_{check["method"]}_c{c}_rep1.strands.txt';mismatch=0
  for seq in path.read_text().splitlines():
   word=dna_to_bytes(seq);seed=int.from_bytes(word[:4],'big')
   if seed not in expected_by_seed:
    _,_,indices=prng.get_src_blocks(seed);value=0
    for index in indices:value^=int.from_bytes(chunks[index],'big')
    expected_by_seed[seed]=value.to_bytes(32,'big')
   mismatch+=word[4:36]!=expected_by_seed[seed]
  check['payload_mismatches_against_reencoded_archive']=mismatch;check['original_synthesis_seed_list_membership_verified']=False
(R/'fountain_file_validation.json').write_text(json.dumps(filechecks,indent=2));(R/'results_validation.json').write_text(json.dumps(dict(completed_runs=27,all_exit_zero=True,repeat_checks=checks,fountain=filechecks),indent=2))
