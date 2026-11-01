"""After timing/audits/plots, author and render the manuscript; visual QA remains a human-visible gate."""
from pathlib import Path
import hashlib,json,subprocess,sys,time
R=Path(__file__).resolve().parent;PROJECT=R.parents[3];INTERNAL=PROJECT/'output/internal';flag=INTERNAL/'resource_docx_operation_started.json'
while True:
 try:
  if 'RESOURCE_AUDITS_AND_FIGURES_COMPLETE' in (R/'postprocess.log').read_text()[-500:]:break
 except FileNotFoundError:pass
 time.sleep(15)
node='/Users/song/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node';skill=Path('/Users/song/.codex/plugins/cache/openai-primary-runtime/documents/26.1015.12040/skills/documents')
if not flag.exists():
 subprocess.run([node,str(skill/'container_tools/mark_artifact_operation_started.mjs'),'--operation-kind','edit','--expected-output-count','1','--output-format','docx'],check=True)
 flag.write_text(json.dumps(dict(operation='edit',expected_output_count=1,format='docx',marker_ran_successfully=True),indent=2))
subprocess.run([sys.executable,str(R/'update_resource_manuscript.py')],check=True)
doc=PROJECT/'output/docx/DBGPS2_manuscript_speed_comparison_serial_3workers.docx';render=Path('/tmp/dbgps2_resource_controls_20261019/render_v4')
subprocess.run([sys.executable,str(skill/'render_docx.py'),str(doc),'--output_dir',str(render),'--emit_pdf'],check=True)
from PIL import Image,ImageDraw
old=Path('/tmp/dbgps2_speed_bloom_20261019/render');oldhash={hashlib.sha256(p.read_bytes()).hexdigest():p.name for p in old.glob('page-*.png')};pages=sorted(render.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[1]));changed=[];inherited=[]
for p in pages:
 match=oldhash.get(hashlib.sha256(p.read_bytes()).hexdigest())
 if match:inherited.append(dict(page=p.name,verified_baseline_page=match))
 else:changed.append(p)
for batch in range(0,len(changed),6):
 canvas=Image.new('RGB',(1530,1320),'#dddddd');draw=ImageDraw.Draw(canvas)
 for i,p in enumerate(changed[batch:batch+6]):
  image=Image.open(p);image.thumbnail((500,625));x=i%3*510;y=i//3*660;canvas.paste(image,(x,y+25));draw.text((x+8,y+5),p.name,fill='black')
 canvas.save(render.parent/f'contact-{batch//6+1}.png')
(INTERNAL/'resource_render_qa_inventory.json').write_text(json.dumps(dict(pages=len(pages),changed_pages_to_inspect=[p.name for p in changed],inherited_pixel_identical_pages=inherited,render=str(render),visual_qa_complete=False),indent=2));print('DOCUMENT_RENDERED_AWAITING_VISUAL_QA',len(pages),len(changed),flush=True)
