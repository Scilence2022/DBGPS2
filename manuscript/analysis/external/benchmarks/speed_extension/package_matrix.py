"""Portable evidence package; no raw reads, QA renders or excluded HT5 results."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED,ZIP_STORED
import hashlib,json
R=Path(__file__).resolve().parent
DEST=R.parent/'DBGPS2_speed_matrix_HT176_Bloom.zip'
def main():
 rows=json.loads((R/'matrix_results.json').read_text());valid=json.loads((R/'matrix_validation.json').read_text());assert len(rows)==81 and valid['completed_runs']==81
 prefix='DBGPS2_speed_matrix_HT176_Bloom/';manifest={};index=[];seen={}
 rootfiles=['LICENSE','requirements.txt','build_matrix_tools.py','docx_ooxml_helpers.py','package_matrix.py','README.md','prepare.py','provenance.py','run_benchmark.py','run_bloom_matrix.py','msa_pipeline.py','smoke_test.py','fragment_smoke_test.py','bloom_smoke_test.py','validate_matrix.py','summarize_matrix.py','update_matrix_manuscript.py','input_provenance.json','input_revalidation.json','annotation_preparation.json','replay_annotation.py','selected_settings.json','matrix_metadata.json','matrix_results.json','matrix_summary.json','matrix_validation.json','fountain_file_validation.json','smoke_validation.json','fragment_smoke_validation.json','bloom_smoke_validation.json','Table6_configuration_comparison.csv','matrix_all_repeats.csv','matrix_stage_timings.csv']
 with ZipFile(DEST,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
  def add(path,name,stored=False):
   path=Path(path);data=path.read_bytes();manifest[name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest());z.writestr(prefix+name,data,compress_type=ZIP_STORED if stored else ZIP_DEFLATED)
  for name in rootfiles:add(R/name,name)
  for pattern in ['Fig10_reconstruction_speed_Bloom_matrix.*','Fig11_Bloom_memory_and_recovery.*']:
   for path in R.glob(pattern):add(path,path.name)
  for dirname in ['src','bin']:
   for path in (R/dirname).rglob('*'):
    if path.is_file() and '__pycache__' not in path.parts:add(path,str(path.relative_to(R)))
  for path in (R/'sources').rglob('*'):
   if path.is_file() and (path.suffix in ['.json','.gz'] or path.name in ['LICENSE','COPYING']):
    if path.name.startswith('figshare_'):continue
    add(path,str(path.relative_to(R)),stored=path.suffix=='.gz')
  for name in ['fountain','song_ht176']:
   for suffix in ['.preparation.json','.merge.stderr']:add(R/'data'/(name+suffix),'preparation/'+name+suffix)
  for row in rows:
   path=Path(row['output_path']);sha=row['sha256']
   if sha not in seen:
    name='canonical_outputs/'+row['dataset']+'_'+sha[:16]+'.txt';add(path,name);seen[sha]=name
   index.append(dict(dataset=row['dataset'],case=row['case'],repeat=row['repeat'],sha256=sha,archive_path=seen[sha],original_path=str(path)))
   result=Path(row['result_path']);stem=str(result)[:-len('.result.json')]
   for suffix in ['.result.json','.stdout','.stderr','.stages.json','.load_before.txt']:
    p=Path(stem+suffix)
    if p.exists():add(p,'runs/'+p.name)
  for path in (R/'runs').glob('fountain_*.decode.log'):add(path,'runs/'+path.name)
  for path in (R/'pilot').glob('*'):
   if path.is_file():add(path,'pilot/'+path.name)
  refs={'song_designs.fa':Path('/Users/song/dbgps2_ht/ref.fa'),'organick_file16_designs.fa':Path('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa')}
  for name,path in refs.items():add(path,'reference_evaluation_only/'+name)
  old=R.parent/'DBGPS2_reconstruction_speed_benchmark.zip'
  if old.exists():add(old,'photolithographic/DBGPS2_reconstruction_speed_benchmark.zip',stored=True)
  z.writestr(prefix+'canonical_output_index.json',json.dumps(index,indent=2));z.writestr(prefix+'MANIFEST.json',json.dumps(manifest,indent=2))
 with ZipFile(DEST) as z:assert z.testzip() is None
 report=dict(path=str(DEST),bytes=DEST.stat().st_size,sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),matrix_entries=81,canonical_unique_outputs=len(seen),excluded=['raw reads','processed FASTA','HT5 measurements','development and rejected runs','document QA renders'])
 (R/'package_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
