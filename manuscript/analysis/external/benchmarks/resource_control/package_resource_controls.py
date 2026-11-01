"""Package completed measurements and reproduction sources without raw inputs or QA renders."""
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib,json
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension';DEST=R.parent/'DBGPS2_serial_and_3worker_controls.zip'
def main():
 rows=json.loads((R/'resource_results.json').read_text());v=json.loads((R/'resource_validation.json').read_text());assert len(rows)==162 and v['completed_runs']==162
 assert json.loads((R/'document_validation.json').read_text())['visual_qa_complete']
 prefix='DBGPS2_resource_controls/';manifest={};seen={};index=[]
 with ZipFile(DEST,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
  def add(p,name):
   data=Path(p).read_bytes();manifest[name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest());z.writestr(prefix+name,data)
  for p in R.glob('*'):
   if p.is_file() and (p.suffix in ['.py','.md','.patch','.json','.csv','.svg','.pdf','.log'] or p.name.startswith('Fig') and p.suffix=='.png') and p.name not in ['package_validation.json','visual_qa.json']:
    add(p,'resource_control/'+p.name)
  for dirname in ['src','bin']:
   for p in (R/dirname).rglob('*'):
    if p.is_file() and '__pycache__' not in p.parts:add(p,'resource_control/'+str(p.relative_to(R)))
  for row in rows:
   p=Path(row['output_path']);sha=row['sha256']
   if sha not in seen:
    name='resource_control/canonical_outputs/'+row['dataset']+'_'+sha[:16]+'.txt';add(p,name);seen[sha]=name
   index.append(dict(workers=row['workers'],dataset=row['dataset'],case=row['case'],repeat=row['repeat'],sha256=sha,archive_path=seen[sha]))
   rp=Path(row['result_path']);stem=str(rp).removesuffix('.result.json')
   for suffix in ['.result.json','.stdout','.stderr','.stages.json','.load_before.txt']:
    p=Path(stem+suffix)
    if p.exists():add(p,'resource_control/'+row['resource_mode']+'/runs/'+p.name)
  # The sibling dependencies preserve relative imports in the new scripts.
  for name in ['run_benchmark.py','run_bloom_matrix.py','docx_ooxml_helpers.py','selected_settings.json','matrix_results.json','input_provenance.json','annotation_preparation.json','fountain_file_validation.json','matrix_validation.json','README.md','LICENSE','requirements.txt']:
   add(B/name,'speed_extension/'+name)
  for p in (B/'src').glob('*'):
   if p.is_file():add(p,'speed_extension/src/'+p.name)
  for name in ['fountain','song_ht176']:add(B/'data'/(name+'.preparation.json'),'speed_extension/data/'+name+'.preparation.json')
  add(R/'bin/muscle3','native_dependencies/bin/muscle3')
  prior=R.parent/'reconstruction_speed'
  for name in ['build_muscle.py','input_and_build_provenance.json','sources/muscle3.8.31_src.tar.gz']:
   add(prior/name,'native_dependencies/'+name)
  for name,path in {'song_designs.fa':Path('/Users/song/dbgps2_ht/ref.fa'),'organick_file16_designs.fa':Path('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa')}.items():add(path,'reference_evaluation_only/'+name)
  z.writestr(prefix+'canonical_output_index.json',json.dumps(index,indent=2));z.writestr(prefix+'MANIFEST.json',json.dumps(manifest,indent=2))
 with ZipFile(DEST) as z:assert z.testzip() is None
 report=dict(path=str(DEST),bytes=DEST.stat().st_size,sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),new_measurements=162,unique_output_sets=len(seen),excluded=['raw FASTQ','processed FASTA','HT5 measurements','development runs','document QA renders'])
 (R/'package_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
