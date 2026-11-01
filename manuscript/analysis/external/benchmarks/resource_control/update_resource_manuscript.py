"""Patch the verified Bloom manuscript after both new controls are complete."""
from pathlib import Path
from copy import deepcopy
from zipfile import ZipFile
import csv,hashlib,importlib.util,json
from lxml import etree as E
from PIL import Image
R=Path(__file__).resolve().parent;B=R.parent/'speed_extension';spec=importlib.util.spec_from_file_location('docxhelpers',B/'docx_ooxml_helpers.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);NS=h.NS;W=h.W
SOURCE=R.parents[2]/'docx/DBGPS2_manuscript_speed_comparison_HT176_Bloom.docx';DEST=SOURCE.with_name('DBGPS2_manuscript_speed_comparison_serial_3workers.docx')
def txt(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def main():
 s=json.loads((R/'resource_summary.json').read_text());validation=json.loads((R/'resource_validation.json').read_text());assert len(s)==54 and validation['completed_runs']==162 and validation['all_outputs_match_original_eight_worker_sets']
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='caaef177bcb3d1788393f921a40cc6b86cc43c2989989b995543c1ad8b10b52a'
 get=lambda w,d,c:next(z for z in s if z['workers']==w and z['dataset']==d and z['case']==c)
 with ZipFile(SOURCE) as z:parts={n:z.read(n) for n in z.namelist()};infos=z.infolist()
 root=E.fromstring(parts['word/document.xml']);body=root.find('w:body',NS);ps=body.findall('w:p',NS);find=lambda prefix:next(p for p in ps if txt(p).startswith(prefix));template=find('We extended the speed comparison to DNA Fountain');heading=find('Reconstruction speed across datasets and Bloom configurations');caption=find('Table 6 Reconstruction');foot=find('Time is mean ± sample s.d. of three complete serial computational repeats')
 old='On the same retained molecules, the tested DBGPS2 configurations reconstructed strands 5.2–14.3-fold faster than this implementation, with three computational repeats per point and format-specific quality endpoints.'
 newer=[z for z in s if z['case'].startswith('dbgps2')];new=f'With fully serial reconstruction or three primary workers, the tested DBGPS2 configurations reconstructed strands {min(z["wall_speedup"] for z in newer):.1f}–{max(z["wall_speedup"] for z in newer):.1f}-fold faster than this implementation on the same retained molecules, with three computational repeats per point and format-specific quality endpoints.'
 h.replace_run_text(find('DNA degradation and sequence errors obstruct data retrieval'),old,new)
 mem=get(3,'song_ht176','dbgps2')['peak_graph_rss_gb'];filtered=min(get(3,'song_ht176',c)['peak_graph_rss_gb'] for c in ['dbgps2_F32','dbgps2_F34'])
 h.replace_run_text(find('DNA degradation and sequence errors obstruct data retrieval'),'Two-pass filtering reduced HT176 peak graph RSS from 12.1 to 6.6 GB, with additional counting work.',f'With three counting workers, two-pass filtering reduced HT176 peak graph RSS from {mem:.1f} to {filtered:.1f} GB at additional counting cost.')
 h.replace_run_text(find('We extended the speed comparison to DNA Fountain'),'We extended the speed comparison','Using eight primary counting or alignment workers, we extended the speed comparison')
 h.append_text(foot,' Table 6 reports the original eight-primary-worker operating points; the serial and three-primary-worker controls are reported separately in Figures 12 and 13 and Table 7.')
 foot_rpr=foot.find('w:r/w:rPr',NS)
 if foot_rpr is not None:foot.findall('w:r',NS)[-1].insert(0,deepcopy(foot_rpr))
 for prefix,sentence in [('Fig. 10 |',' Counting and alignment use eight primary workers in this figure.'),('Fig. 11 |',' This figure uses eight graph insertion workers with the original pipeline.')]:
  cp=find(prefix);h.append_text(cp,sentence);rp=E.Element(W+'rPr');rp.append(h.element('sz',val=22));cp.findall('w:r',NS)[-1].insert(0,rp)
 nodes=[h.paragraph('Reconstruction with one serial worker and three primary workers',heading)]
 nodes.append(h.paragraph('We repeated the complete three-dataset, nine-configuration matrix under two additional resource protocols to distinguish operating-system process count from computational parallelism. In the serial protocol, one pipeline worker performed read parsing, k-mer extraction and hash insertion sequentially, with counting callbacks executed inline on that worker; seed enumeration, filter construction and path search were also serial. MSA processed one cluster at a time without a Python process pool, waiting for its single-threaded MUSCLE child before voting and starting the next cluster. This protocol limits algorithmic parallelism to one active reconstruction worker, without pinning that worker to one physical core.',template))
 nodes.append(h.paragraph('The three-primary-worker protocol used three graph insertion workers with the existing three-stage reading/extraction/counting pipeline, and three MSA workers with a coordinator. Each MSA worker waited for its single-threaded MUSCLE child and performed voting between calls. These are matched primary worker counts, not a strict limit of three total threads or three CPU cores. All inputs, search parameters, read selection, validity checks and output rules were retained. The 162 new measurements comprise three complete computational repeats of nine configurations on three inputs under each protocol; no earlier timing was reused (Fig. 12).',template))
 for w,word in [(1,'fully serial'),(3,'three-primary-worker')]:
  db2=[z for z in s if z['workers']==w and z['case'].startswith('dbgps2')];primary=[get(w,d,'dbgps2') for d in ['fountain','nbt_file16','song_ht176']]
  nodes.append(h.paragraph(f'Under the {word} protocol, wall-time speed ratios for the tested DBGPS2 configurations ranged from {min(z["wall_speedup"] for z in db2):.1f} to {max(z["wall_speedup"] for z in db2):.1f} relative to the corresponding clustering/MSA/voting measurement. Primary no-Bloom DBGPS2 speed ratios were '+', '.join(f'{z["wall_speedup"]:.1f}' for z in primary)+' for DNA Fountain, Organick File 16 and HT176, respectively. Table 7 reports wall and user-plus-system CPU times for the primary methods; the complete configuration matrix and its quality endpoints accompany the source data.',template))
 nodes.append(h.paragraph('All 162 completed outputs matched the corresponding previously verified eight-worker output sets by SHA256, including agreement across resource settings and computational repeats. Query filtering again preserved the primary no-Bloom outputs; F32, F34 and B12/F34 agreed with one another. The existing output counts and fidelity endpoints in Table 6 therefore also apply to these scheduling controls. For Fountain, identical output hashes link the new timings to the independently checksum-verified archive and payload re-encoding checks. Timing comparisons retain their implementation-specific and operating-point scope.',template))
 nodes.append(h.paragraph('Table 7 Reconstruction wall and CPU time under serial and three-primary-worker protocols',caption))
 rows=[['Input','Workers','Method','Wall time (s)','CPU time (s)','Speedup']]
 for d,label in [('fountain','DNA Fountain'),('nbt_file16','Organick File 16'),('song_ht176','HT176')]:
  for w in [1,3]:
   for case,name in [('cluster_msa_vote','Clustering + MSA + vote'),('dbgps_rule','DBGPS rule (adapted)'),('dbgps2','DBGPS2 no Bloom')]:
    z=get(w,d,case);rows.append([label,str(w),name,f'{z["wall_mean"]:,.1f} ± {z["wall_sd"]:.1f}',f'{z["cpu_mean"]:,.1f} ± {z["cpu_sd"]:.1f}',f'{z["wall_speedup"]:.1f}×'])
 table=h.table(rows);widths=[1200,1050,1750,2150,2150,1060]
 for col,width in zip(table.findall('w:tblGrid/w:gridCol',NS),widths):col.set(W+'w',str(width))
 for row in table.findall('w:tr',NS):
  for ci,(cell,width) in enumerate(zip(row.findall('w:tc',NS),widths)):
   cell.find('w:tcPr/w:tcW',NS).set(W+'w',str(width));cell.find('w:p/w:pPr/w:jc',NS).set(W+'val','left' if ci==2 else 'center')
 nodes.append(table);nodes.append(h.paragraph('Values are mean ± sample s.d. across three newly measured computational repeats. CPU time sums user and system time of reconstruction and waited-for children, plus calling-harness self CPU through final output canonicalization. Wall time includes process startup, input parsing, reconstruction, validity filtering, conflict rejection and output writing. Shared preprocessing, parameter tuning, reference scoring, hash comparisons and outer decoding are excluded. Speedup is the ratio of the corresponding mean MSA-baseline wall time to mean method wall time within the same input and resource protocol. Primary search parameters and all quality endpoints are unchanged from Table 6. The three-worker protocol matches insertion/alignment workers and retains pipeline and coordination overhead; it is not a strict three-core allocation.',foot))
 # Describe phase measurements with full numerical source tables supplied separately.
 stages=list(csv.DictReader((R/'stage_wall_and_cpu.csv').open()));phase=lambda w,d,c,p:next(float(z['wall_mean']) for z in stages if int(z['workers'])==w and z['dataset']==d and z['case']==c and z['stage']==p)
 for d,label in [('fountain','DNA Fountain'),('nbt_file16','Organick File 16'),('song_ht176','HT176')]:
  nodes.append(h.paragraph(f'For {label}, the MSA/voting stage took '+ ' and '.join(f'{phase(w,d,"cluster_msa_vote","msa_vote"):.1f} s' for w in [1,3])+' under the serial and three-primary-worker protocols, respectively. Primary DBGPS2 counting took '+ ' and '.join(f'{phase(w,d,"dbgps2","count"):.1f} s' for w in [1,3])+', seed handling '+ ' and '.join(f'{phase(w,d,"dbgps2","seed"):.1f} s' for w in [1,3])+', and search with validity checking '+ ' and '.join(f'{phase(w,d,"dbgps2","search"):.1f} s' for w in [1,3])+'. Source tables report the wall and CPU time of every stage for every configuration, including both counting passes and query-filter construction.',template))
 nodes.append(h.paragraph('Figure 13 distinguishes elapsed time from CPU consumption. Multiple workers can reduce elapsed time without an equivalent reduction in total CPU time; process creation, communication and synchronization contribute to the latter. Computational repeat variability and the retained background-process records describe these workstation measurements. The serial control removes concurrent reconstruction workers, while the three-worker control evaluates the requested primary worker setting. Neither protocol fixes affinity, equal physical-core placement or an otherwise idle host.',template))
 ix=body.index(find('Discussion'))
 for n in nodes:body.insert(ix,n);ix+=1
 h.append_text(find('The expanded comparison places reconstruction mechanisms'),' The newly measured serial and three-primary-worker controls (Figs. 12 and 13) additionally distinguish worker scheduling from process count and report CPU consumption together with elapsed time. Their output sets were identical to the original eight-worker operating points; CPU and stage measurements accompany each new timing rather than being inferred by multiplying elapsed time by a configured worker count.')
 methods=[
 'For the resource controls, benchmark-only scheduling changes exposed graph count-worker and pipeline-worker counts without changing the counter, Bloom operations, traversal or validators. The serial path used kt_pipeline(1) and kt_for(1), with non-overlapping read/extract/insert stages; the latter executes its callback inline. The three-primary-worker path used kt_pipeline(3) and kt_for(3). Search, seed handling and query-filter construction remained serial. MSA used a direct sequential map for one worker and the original ProcessPoolExecutor with three workers for the other protocol. Numerical-library parallelism was limited to one. Fifty format, Bloom, serial/parallel and fragmented-read smoke checks passed before real-data timings.',
 'All resource-control measurements were new, complete runs scheduled serially with cyclic configuration order and alternating resource-setting order within input/repeat blocks. Eight-worker output hashes were consulted only after each timer ended to test scheduling invariance. Raw inputs and references were independently rehashed after all measurements. No forced cold-cache regime, CPU affinity or dedicated idle host was used; background-process snapshots were retained. The three-worker setting matches primary insertion and alignment workers but also includes graph read/extract stages and MSA process-pool coordination, so it does not enforce equal total thread counts or a strict three-core budget.',
 'Complete reconstruction wall time spans native process startup through validity checks, reference-free conflict rejection and sorted output writing. User and system CPU time from /usr/bin/time -l includes reconstruction and waited-for children; calling-harness self CPU over final output canonicalization was added and retained separately. Native graph stage CPU time was sampled with getrusage(RUSAGE_SELF), which includes its threads. MSA phase CPU combines self and reaped-child usage, sampled after the worker pool has closed and joined. Reported phases are count, query-filter construction, seed and search for graphs, or clustering, MSA/voting and validation for the conventional pipeline. Process startup/shutdown, graph destruction and final harness work contribute to the complete timer but not all phase timers; source tables retain the residual explicitly. CPU time is measured consumption rather than an estimate based on requested workers.'
 ]
 ix=body.index(find('All 81 matrix entries are complete measurements'))+1
 for text in methods:body.insert(ix,h.paragraph(text,template));ix+=1
 h.append_text(find('Figure 8 and Supplementary Figs. S3–S5 were generated'),' Figures 12 and 13 and Table 7 use the 162 newly measured resource-control runs; source CSVs report every repeat, configuration, wall/CPU phase, output hash and resource setting.')
 h.append_text(find('DBGPS2 is available at https://github.com/'),' The resource_control package supplies the serial and three-worker scheduling adapters, instrumentation, 162-run records, correctness checks and source tables and vector figures for the new controls.')
 rels=E.fromstring(parts['word/_rels/document.xml.rels']);pic_template=next(p for p in reversed(ps) if p.find('.//w:drawing',NS) is not None);legend_template=find('Fig. 11 |');images={}
 def figure(num,name,legend):
  path=R/name;pic=deepcopy(pic_template);ppr=pic.find('w:pPr',NS)
  for tag in ['pageBreakBefore','keepNext']:
   old=ppr.find('w:'+tag,NS)
   if old is not None:ppr.remove(old)
   ppr.append(h.element(tag))
  im=Image.open(path);cx=5943600;cy=int(cx*im.height/im.width)
  for el in [pic.find('.//wp:extent',NS),*pic.findall('.//a:ext',NS)]:el.set('cx',str(cx));el.set('cy',str(cy))
  prop=pic.find('.//wp:docPr',NS);prop.set('id',str(max(int(x.get('id')) for x in root.findall('.//wp:docPr',NS))+1));prop.set('name',f'Figure {num}');prop.set('descr',legend);pic.find('.//pic:cNvPr',NS).set('name',name)
  rid='rId'+str(max(int(x.get('Id')[3:]) for x in rels if x.get('Id').startswith('rId'))+1);E.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',Id=rid,Type=NS['r']+'/image',Target='media/'+name)
  pic.find('.//a:blip',NS).set('{'+NS['r']+'}embed',rid)
  p=h.paragraph(legend,legend_template);pp=p.find('w:pPr',NS)
  for tag in ['keepNext','pageBreakBefore']:
   old=pp.find('w:'+tag,NS)
   if old is not None:pp.remove(old)
   pp.append(h.element(tag,val=0))
  split=legend.index('. ',5);lead,tail=legend[:split+1],legend[split+2:];h.replace_text(p,lead);h.append_text(p,' '+tail)
  for i,run in enumerate(p.findall('w:r',NS)):
   rp=E.Element(W+'rPr');rp.append(h.element('sz',val=22))
   if i==0:rp.append(h.element('b'))
   run.insert(0,rp)
  ix=len(body)-1 if body[-1].tag==W+'sectPr' else len(body);body.insert(ix,pic);body.insert(ix+1,p);images['word/media/'+name]=path.read_bytes()
 figure(12,'Fig12_resource_control_speed.png','Fig. 12 | Reconstruction speed under serial and three-primary-worker protocols. Rows show DNA Fountain, Organick File 16 and HT176; the left column uses serial reconstruction and the right column uses three primary workers. Each panel compares clustering/MSA/voting, the adapted DBGPS rule and seven DBGPS2 configurations. Bars show mean complete reconstruction wall time with sample s.d. across three new computational repeats. Axes are logarithmic; annotations are wall-time ratios to the MSA baseline in the same input/protocol. B8/B16 denote query filters; F32/F34 denote two-pass counting filters; B12/F34 combines both. Serial graph counting uses one pipeline worker with inline insertion; MSA handles one cluster at a time. The three-worker graph protocol retains reading/extraction stages and the MSA protocol retains its coordinator, so neither is described as a strict three-core limit. All search parameters and output sets match Table 6. Timers include filter construction or both counting passes, validity checks and final output handling; common preparation, scoring, hash comparison and outer decoding are excluded.')
 figure(13,'Fig13_resource_control_wall_cpu.png','Fig. 13 | Complete reconstruction wall time and CPU consumption. Primary clustering/MSA/voting (MSA), adapted DBGPS (DBG) and no-Bloom DBGPS2 (DB2) are shown for serial (1) and three-primary-worker (3) protocols on the three inputs. Blue bars show elapsed wall time and orange bars show measured user-plus-system CPU time, including waited-for children and final harness output handling. Bars are means with sample s.d. across three new computational repeats; axes are logarithmic. Table 7 gives the numeric values. CPU time was measured rather than inferred from worker count. Stage wall and CPU times for every method, Bloom configuration, input and repeat accompany the source data. Serial reconstruction has one active algorithmic worker; the three-worker setting retains graph pipeline and process-coordination overhead and does not impose a strict physical-core allocation.')
 # Distinguish sequential repeat scheduling from the new serial resource protocol.
 for p in root.findall('.//w:p',NS):
  if 'serial computational repeats' in txt(p):h.replace_run_text(p,'serial computational repeats','sequential computational repeats')
 updates={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),'word/_rels/document.xml.rels':E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)}
 with ZipFile(DEST,'w') as z:
  for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
  for name,data in images.items():z.writestr(name,data)
 with ZipFile(DEST) as z:
  assert z.testzip() is None
  for name,data in parts.items():
   if name.startswith('word/media/'):assert z.read(name)==data
 oldroot=E.fromstring(parts['word/document.xml']);m={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'};assert [E.tostring(x) for x in oldroot.findall('.//m:oMath',m)]==[E.tostring(x) for x in root.findall('.//m:oMath',m)]
 report=dict(output=str(DEST),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),figures=len(root.findall('.//w:drawing',NS)),tables=len(root.findall('.//w:tbl',NS)),new_measurements=162,preserved_old_media_and_equations=True)
 (R/'document_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
