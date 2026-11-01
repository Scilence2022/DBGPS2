"""Run heavy audits only after the entire measured matrix has completed."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time
R=Path(__file__).resolve().parent
while True:
 try:
  rows=json.loads((R/'matrix_results.json').read_text())
  if len(rows)==81:break
 except (FileNotFoundError,json.JSONDecodeError):pass
 time.sleep(3)
time.sleep(1)
print('All 81 measurements completed; starting untimed audits',flush=True)
inputs=json.loads((R/'input_provenance.json').read_text());checks=[]
for row in inputs:
 p=Path(row['path']);sh=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8*1024*1024):sh.update(b)
 assert p.stat().st_size==row['bytes'] and sh.hexdigest()==row['sha256'],str(p)
 checks.append(dict(path=str(p),sha256=sh.hexdigest(),unchanged=True))
(R/'input_revalidation.json').write_text(json.dumps(checks,indent=2));print('Input hashes unchanged',flush=True)
# Record the graph-only synthetic boundary-view preparation separately.
# This replay occurs after all reconstruction timers; no payload bases change.
annotation_tmp=R/'data/fountain.graph.annotation_replay.fa'
annotation_start=time.perf_counter()
with (R/'data/fountain.fa').open('rb') as f,annotation_tmp.open('wb') as out:
 while header:=f.readline():
  seq=f.readline().strip();out.write(header+b'CCTGCAGAGT'+seq+b'T\n')
annotation_seconds=time.perf_counter()-annotation_start
sh=hashlib.sha256()
with annotation_tmp.open('rb') as f:
 while b:=f.read(8*1024*1024):sh.update(b)
expected=next(x['sha256'] for x in inputs if x['path'].endswith('/fountain.graph.fa'))
assert sh.hexdigest()==expected
annotation_tmp.unlink()
(R/'annotation_preparation.json').write_text(json.dumps(dict(dataset='fountain',synthetic_prefix='CCTGCAGAGT',synthetic_terminal_lookahead='T',wall_seconds=annotation_seconds,repeats=1,measured_after_reconstruction_matrix=True,sha256_matches_original_graph_view=True,payload_sequence_unchanged=True,in_reconstruction_timer=False),indent=2))
subprocess.run([sys.executable,str(R/'validate_matrix.py')],check=True)
env=os.environ.copy();env['PYTHONPATH']='/tmp/dbgps2_speed_plot_dependencies:/tmp/dbgps2_speed_extension_dependencies'
subprocess.run([sys.executable,str(R/'summarize_matrix.py')],check=True,env=env)
print('AUDITS_AND_FIGURES_COMPLETE',flush=True)
