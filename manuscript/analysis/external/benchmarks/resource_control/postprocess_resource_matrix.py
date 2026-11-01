"""Wait without competing for resources, then audit inputs and produce source tables."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,os
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension'
while True:
 try:
  rows=json.loads((R/'resource_results.json').read_text());log=(R/'resource_matrix.log').read_text()
  if len(rows)==162 and 'RESOURCE_MATRIX_ALL_DONE 162' in log[-300:]:break
 except (FileNotFoundError,json.JSONDecodeError):pass
 time.sleep(15)
print('TIMED_MATRIX_ENDED_STARTING_AUDITS',flush=True)
checks=[]
for old in json.loads((B/'input_provenance.json').read_text()):
 p=Path(old['path']);h=hashlib.sha256()
 with p.open('rb') as f:
  while block:=f.read(8*1024*1024):h.update(block)
 assert h.hexdigest()==old['sha256'] and p.stat().st_size==old['bytes'],p
 checks.append(dict(path=str(p),sha256=h.hexdigest(),unchanged=True));print('INPUT_UNCHANGED',p.name,flush=True)
(R/'input_revalidation.json').write_text(json.dumps(checks,indent=2))
manifest=json.loads((R/'source_manifest.json').read_text())
for name,old in manifest.items():
 p=R/name;assert hashlib.sha256(p.read_bytes()).hexdigest()==old['sha256'],p
print('SOURCES_AND_BINARIES_UNCHANGED',flush=True)
subprocess.run([sys.executable,str(R/'validate_and_summarize.py')],check=True)
if (R/'plot_resource_matrix.py').exists():
 env=os.environ.copy();env['PYTHONPATH']='/tmp/dbgps2_speed_plot_dependencies:/tmp/dbgps2_speed_extension_dependencies';subprocess.run([sys.executable,str(R/'plot_resource_matrix.py')],env=env,check=True)
print('RESOURCE_AUDITS_AND_FIGURES_COMPLETE',flush=True)
