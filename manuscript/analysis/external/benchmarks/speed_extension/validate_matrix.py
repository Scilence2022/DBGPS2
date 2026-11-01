"""Audit all completed timings, Bloom output invariants and untimed file recovery."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
from run_bloom_matrix import CASES
R=Path(__file__).resolve().parent
rows=json.loads((R/'matrix_results.json').read_text());assert len(rows)==81
checks=[]
for dataset in ['fountain','nbt_file16','song_ht176']:
 for case in CASES:
  group=[z for z in rows if z['dataset']==dataset and z['case']==case];assert len(group)==3 and {z['repeat'] for z in group}=={1,2,3};assert len({z['sha256'] for z in group})==1;assert len({(z['emitted'],z['exact'],z['incorrect']) for z in group})==1
  for z in group:
   p=Path(z['output_path']);b=p.read_bytes();assert z['exit_code']==0 and z['wall_seconds']>0 and hashlib.sha256(b).hexdigest()==z['sha256'] and len(b.splitlines())==z['emitted']
  checks.append(dict(dataset=dataset,case=case,repeats=3,identical=True,sha256=group[0]['sha256']))
 get=lambda c:next(z for z in rows if z['dataset']==dataset and z['case']==c)
 for c in ['dbgps2_B8','dbgps2_B16']:assert get(c)['sha256']==get('dbgps2')['sha256'],(dataset,c,'query filter changed output')
 for c in ['dbgps2_F34','dbgps2_B12_F34']:assert get(c)['sha256']==get('dbgps2_F32')['sha256'],(dataset,c,'counting filter sets differ')
 print('VERIFIED',dataset,'nine configurations, three repeats',flush=True)
# Full-file decoding is outside every reconstruction timer; cache identical sets.
expected='8651e90d3a013178b816b63fdbb94b9b';decoded_by_hash={};filechecks=[]
for case in CASES:
 row=next(z for z in rows if z['dataset']=='fountain' and z['case']==case and z['repeat']==1);sha=row['sha256']
 if sha not in decoded_by_hash:
  out=R/'runs'/f'fountain_{case}.decoded.bin';cmd=[sys.executable,str(R/'src/decode_fountain.py'),row['output_path'],'--no_screen','--out',str(out)];env=os.environ.copy();env['PYTHONPATH']='/tmp/dbgps2_speed_extension_dependencies';cp=subprocess.run(cmd,capture_output=True,text=True,env=env);(R/'runs'/f'fountain_{case}.decode.log').write_text(cp.stdout+cp.stderr);md5=hashlib.md5(out.read_bytes()).hexdigest() if out.exists() else None
  decoded_by_hash[sha]=dict(command=cmd,exit_code=cp.returncode,bytes=out.stat().st_size if out.exists() else None,md5=md5,md5_match=md5==expected,decoded_path=str(out),expected_md5=expected)
 check=dict(case=case,method=case,output_sha256=sha,decoding_in_timing=False,**decoded_by_hash[sha]);filechecks.append(check);print('FOUNTAIN_FILE',case,check['md5_match'],flush=True)
# Independently re-encode every accepted seed/payload against a verified archive.
verified=next((x for x in filechecks if x['md5_match']),None)
if verified:
 sys.path.insert(0,str(R/'src'));sys.path.insert(0,'/tmp/dbgps2_speed_extension_dependencies')
 from fountain_common import PRNG,dna_to_bytes,ReedSolomon2
 blob=Path(verified['decoded_path']).read_bytes();chunks=[int.from_bytes(blob[i:i+32],'big') for i in range(0,len(blob),32)];prng=PRNG(67088,.001,.025);rs=ReedSolomon2();seed_cache={};word_cache={}
 for check in filechecks:
  row=next(z for z in rows if z['dataset']=='fountain' and z['case']==check['case'] and z['repeat']==1);mismatch=invalid=0
  for seq in Path(row['output_path']).read_text().splitlines():
   if seq not in word_cache:
    word=dna_to_bytes(seq);seed=int.from_bytes(word[:4],'big')
    if seed not in seed_cache:
     prng.set_seed(seed);_,_,indices=prng.get_src_blocks();v=0
     for index in indices:v^=chunks[index]
     seed_cache[seed]=v.to_bytes(32,'big')
    word_cache[seq]=(word[4:36]!=seed_cache[seed],bytes(rs.encode(word[:36]))!=word)
   a,b=word_cache[seq];mismatch+=a;invalid+=b
  assert invalid==0,(check['case'],invalid)
  check.update(payload_mismatches_against_reencoded_archive=mismatch,independent_rs_reencoding_mismatches=invalid,original_synthesis_seed_list_membership_verified=False)
  print('FOUNTAIN_PAYLOAD',check['case'],mismatch,invalid,flush=True)
(R/'fountain_file_validation.json').write_text(json.dumps(filechecks,indent=2));(R/'matrix_validation.json').write_text(json.dumps(dict(completed_runs=81,all_exit_zero=True,repeat_checks=checks,query_filters_preserve_outputs=True,counting_filters_32_34_combined_agree=True,fountain=filechecks),indent=2))
