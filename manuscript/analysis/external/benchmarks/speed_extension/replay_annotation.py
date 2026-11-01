from pathlib import Path
import hashlib,json,time
R=Path(__file__).resolve().parent
inputs=json.loads((R/'input_provenance.json').read_text())
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
