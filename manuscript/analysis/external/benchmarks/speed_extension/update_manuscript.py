"""Extend the verified manuscript in place within its OOXML package."""
from copy import deepcopy
from pathlib import Path
from zipfile import ZipFile
import hashlib,importlib.util,json
from lxml import etree as E
from PIL import Image
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('previous_patch',R.parent/'reconstruction_speed/update_manuscript.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
NS=h.NS;W=h.W
SOURCE=R.parents[2]/'docx/DBGPS2_manuscript_speed_comparison.docx'
DEST=SOURCE.with_name('DBGPS2_manuscript_speed_comparison_expanded.docx')
def txt(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def replace(p,s):h.replace_text(p,s)
def main():
 prep=json.loads((R/'data/song_ht176.preparation.json').read_text());selected=json.loads((R/'selected_settings.json').read_text())
 summary=json.loads((R/'summary.json').read_text());g=lambda d,m:next(z for z in summary if z['dataset']==d and z['method']==m)
 f=lambda z:f'{z["time_mean"]:.1f} ± {z["time_sd"]:.1f}'
 with ZipFile(SOURCE) as z:parts={n:z.read(n) for n in z.namelist()};infos=z.infolist()
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='b3125b15c82b5d71de247050ade1d489ad3922766102fdf72c066c43a4c5c081'
 root=E.fromstring(parts['word/document.xml']);body=root.find('w:body',NS);ps=body.findall('w:p',NS)
 find=lambda s:next(p for p in ps if txt(p).startswith(s))
 intro=find('We evaluate DBGPS2 in two distinct settings.')
 old='Controlled repeated timings additionally compare both graph pruning strategies with prefix clustering, multiple alignment and weighted majority voting on photolithographic Files 1 and 2.'
 h.replace_run_text(intro,old,'Controlled repeated timings compare both graph pruning strategies with prefix clustering, multiple alignment and weighted majority voting on two photolithographic inputs, DNA Fountain, the large Organick File 16 library and the 176-day thermal-aging sample using the Song strand design. The timing and recovery endpoints are stated separately for each format.')
 heading=find('Controlled speed comparison with clustering and multiple alignment');template=find('We compared both graph pruning strategies with the prefix-clustering')
 nodes=[h.paragraph('Speed comparison across additional DNA storage datasets',heading)]
 nodes.append(h.paragraph(f'We extended the benchmark to DNA Fountain (72,000 synthesized droplets), Organick File 16 (607,150 designed strands; nbt.4079), and a 176-day aging sample (HT176) of the 210,000-strand pool using the Song encoding design (Fig. 10 and Table 6). HT176 was aged at 70°C and sequenced in a separate run; local raw-file checksums define its input provenance. The sample uses the published strand-design format; the present timing measurement is specific to HT176. All 15,719,903 raw Organick reads were retained. Shared NGmerge preprocessing yielded 10,477,842 length-152 DNA Fountain molecules and {prep["retained_records"]:,} HT176 fragments; unsuccessful merges were excluded equally for all methods. The HT176 timing experiment therefore uses a common merged-input subset.',template))
 for d,label in [('fountain','DNA Fountain'),('nbt_file16','Organick File 16'),('song_ht176','Song HT176')]:
  a,b,c=[g(d,m) for m in ['cluster_msa_vote','dbgps_rule','dbgps2']]
  s=f'For {label}, clustering/MSA/voting required {f(a)} s, the adapted DBGPS rule {f(b)} s and DBGPS2 {f(c)} s (mean ± s.d., three serial computational repeats). Relative to this conventional implementation, the two graph modes were {b["speedup"]:.1f}- and {c["speedup"]:.1f}-fold faster, respectively. '
  if d=='fountain':
   s+=f'The corresponding counts were {a["emitted"]:,}, {b["emitted"]:,} and {c["emitted"]:,} RS-zero-syndrome droplets with distinct seeds. These are code-valid counts rather than an exact comparison against the original synthesized sequence list. '
   q=json.loads((R/'fountain_file_validation.json').read_text());ok=[x['method'] for x in q if x['md5_match']]
   if len(ok)==3:s+='Each method also recovered the complete archive with the expected MD5 in a separate untimed outer-decoding check.'
   else:s+='Separate outer-decoding outcomes are retained in the source data; complete-file recovery is not inferred from code-valid counts.'
  else:s+=f'Exact encoded-region recovery was {100*a["recovery"]:.2f}%, {100*b["recovery"]:.2f}% and {100*c["recovery"]:.2f}%, with {a["incorrect"]:,}, {b["incorrect"]:,} and {c["incorrect"]:,} incorrect outputs, respectively. '
  nodes.append(h.paragraph(s,template))
 # Reuse the manuscript's numbered bibliography without introducing new entries.
 for marker,reference in [('72,000 synthesized droplets)','4'),('nbt.4079)','7'),('Song encoding design','17')]:
  p=nodes[1]
  run=next(r for r in p.findall('w:r',NS) if marker in ''.join(r.xpath('.//w:t/text()',namespaces=NS)))
  value=''.join(run.xpath('.//w:t/text()',namespaces=NS));before,after=value.split(marker,1)
  run.find('w:t',NS).text=before+marker
  cr=E.Element(W+'r');rp=E.SubElement(cr,W+'rPr');rp.append(h.element('vertAlign',val='superscript'));E.SubElement(cr,W+'t').text=reference
  tail=deepcopy(run);tail.find('w:t',NS).text=after;tail.find('w:t',NS).set('{http://www.w3.org/XML/1998/namespace}space','preserve')
  ix=p.index(run)+1;p.insert(ix,cr);p.insert(ix+1,tail)
 h.replace_run_text(find('For the new comparison, all methods consumed'), 'For the new comparison', 'For the photolithographic comparison')
 caption=h.paragraph('Table 6 Reconstruction time and strand outputs on additional DNA storage datasets',find('Table 5 Reconstruction time'))
 rows=[['Input','Method','Time (s)','Exact / valid','Wrong','Speedup']]
 for d,short in [('fountain','DNA Fountain'),('nbt_file16','Organick File 16'),('song_ht176','Song HT176')]:
  for m,label in [('cluster_msa_vote','Clustering + MSA + vote'),('dbgps_rule','DBGPS rule (adapted)'),('dbgps2','DBGPS2 (adapted)')]:
   z=g(d,m);quality=f'{z["emitted"]:,} RS valid' if d=='fountain' else f'{100*z["recovery"]:.2f}%';wrong='n.a.' if z['incorrect'] is None else str(z['incorrect'])
   rows.append([short,label,f(z),quality,wrong,f'{z["speedup"]:.1f}×'])
 table=h.table(rows);widths=[1300,2000,1800,1700,1000,1560]
 for col,width in zip(table.findall('w:tblGrid/w:gridCol',NS),widths):col.set(W+'w',str(width))
 for row in table.findall('w:tr',NS):
  for cell,width in zip(row.findall('w:tc',NS),widths):cell.find('w:tcPr/w:tcW',NS).set(W+'w',str(width))
 foot=h.paragraph(f'Times exclude shared preprocessing, tuning, reference scoring and outer decoding; the common merge-and-conversion stages took 46.3 s for DNA Fountain and {prep["preparation_seconds"]:.1f} s for HT176 in separate single runs. These costs would be added to each method when budgeting from raw paired FASTQ. Exact recovery uses all designs as denominator for Organick and HT176. Fountain reports RS-valid unique seeds; n.a. indicates that exact erroneous-output counts cannot be established without the original design list. Speedup is the ratio of mean reconstruction times. New-format MSA uses one alignment of up to 15 sampled reads per eligible cluster, whereas Table 5 retains its original five-subset protocol. Graph modes share k within a dataset and have Bloom filters disabled. Outputs and quality counts were identical across repeats.',find('Time is mean ± sample s.d.'))
 nodes.extend([caption,table,foot]);index=body.index(find('Discussion'))
 for n in nodes:body.insert(index,n);index+=1
 disc=find('The controlled reconstruction benchmark extends the resource analysis')
 replace(disc,'The controlled reconstruction benchmarks extend the resource analysis across five processed inputs and several encoding schemes. Avoiding an independent multiple alignment for each prefix cluster can reduce reconstruction cost, but the size of the benefit depends on the input, implementation and selected recovery endpoint (Figs. 9 and 10). This inherited graph-workflow advantage should be distinguished from the incremental effect of the ratio constraint. The two graph modes can trade additional traversal time for improved recovery when their cutoffs are selected separately. The new HT176 benchmark contains only successfully merged molecules, so its time and recovery cannot be substituted for the full-read thermal-pool measurements. Prefix clustering and the per-cluster MUSCLE subprocesses also differ from similarity clustering, in-process alignment and the original Organick BMA decoder. The results establish operating-point comparisons for these implementations; they do not establish a universal speedup over all clustering methods, equal-recovery decoding or complete outer decoding.')
 # Format-specific details appended to timing Methods.
 methodanchor=find('Three repeats per method and input were run sequentially')
 methodnodes=[h.paragraph('The three additional datasets were timed using the same serial repeat and order-rotation scheme. DNA Fountain and HT176 paired reads were merged with NGmerge 0.3 (eight threads, dovetail detection, default overlap and mismatch thresholds). Only merged molecules were retained; Fountain additionally required 152 nt. Shared merging and FASTQ-to-FASTA conversion were measured once and excluded from every reconstruction timer. Organick used the full untrimmed read FASTA. Graph input for Fountain added a constant primer anchor and a one-base boundary lookahead for the existing search interface; these are format metadata and do not supply payload bases. Other methods received the same molecular sequences without these interface constants.',template),
 h.paragraph(f'Observed Fountain 16-nt seed prefixes and forward/reverse Organick primer-plus-prefix seeds were discovered without references. HT176 used the legacy two-base mapping and its fixed numeric index range, without a list of true addresses. Both graph modes used k = 21, 29 and 31 for Fountain, Organick and HT176, respectively. DBGPS2 used (c, r) = (0, 6), (0, 1.9) and (1, 2.55). The absolute-cutoff comparator used r = 1,000,000; c was selected from the fixed grid 1, 3, 5, 10 and 20 (plus 25 and 40 for HT176) on a deterministic hashed 1/128 subset of observed seeds (or enumerated HT176 indices), maximizing exact output count for Organick/HT176 and RS-valid count for Fountain. Selected cutoffs were 1 for Fountain, 3 for Organick and {selected["song_ht176"]["dbgps_cutoff"]} for HT176. This retrospective pilot is excluded from timing and provides no held-out estimate of parameter generalization.',template),
 h.paragraph('The additional-format baseline clustered by exact observed prefixes: 16 nt for Fountain and HT176 and 14 nt after normalizing Organick reads to one orientation. It required three observations, restricted HT176 prefix groups to the same numeric index range as graph search and retained up to 15 reads (seed 20261019), performed one MUSCLE 3.8.31 alignment per cluster, and used gap weight 0.4 in majority voting. Fountain and Organick used uniform reservoir sampling. For fragmented HT176 reads, the baseline preferred reads with the longest observed coverage of the known 164-nt region, used random tie-breaking, and treated terminal alignment gaps as missing observations rather than deletions; internal gaps retained weight 0.4. Primer orientation and trimming used sequence-format metadata. Fountain candidates required zero-syndrome RS(38,36), Organick candidates the same rotating-ternary predicate as graph candidates, and HT176 candidates the same legacy mapping and CRC16. Conflicting candidates for a prefix were rejected without consulting references. Timers included format-specific clustering or counting/seeding/search, validation, conflict rejection and sorted output writing. Exact reference scoring and Fountain outer decoding followed the timer. All repeated output hashes were checked for equality.',template)]
 index=body.index(methodanchor)+1
 for n in methodnodes:body.insert(index,n);index+=1
 p=find('Figure 8 and Supplementary Figs. S3–S5 were generated')
 h.append_text(p,' Figure 10 and Table 6 extend the serial benchmark to three further published libraries; the Fountain panel uses RS-valid seed counts rather than exact design recovery.')
 p=find('DBGPS2 is available at https://github.com/')
 h.append_text(p,' The speed_extension package provides additional input checksums, NGmerge provenance, format adapters, serial repeat logs, output hashes, full-file validation records and source data for Figure 10 and Table 6.')
 # Reuse the existing picture paragraph structure and add one independent part.
 image=R/'Fig10_extended_reconstruction_speed.png';picture=deepcopy(next(p for p in reversed(ps) if p.find('.//w:drawing',NS) is not None))
 for child in list(picture):
  if child.tag not in (W+'pPr',W+'r'):picture.remove(child)
 ppr=picture.find('w:pPr',NS)
 if ppr is None:ppr=E.Element(W+'pPr');picture.insert(0,ppr)
 for tag in ['pageBreakBefore','keepNext']:
  if ppr.find('w:'+tag,NS) is None:ppr.append(h.element(tag))
 width=5943600;height=int(width*Image.open(image).height/Image.open(image).width)
 extent=picture.find('.//wp:extent',NS);extent.set('cx',str(width));extent.set('cy',str(height))
 for ext in picture.findall('.//a:ext',NS):ext.set('cx',str(width));ext.set('cy',str(height))
 dp=picture.find('.//wp:docPr',NS);dp.set('id',str(max(int(x.get('id')) for x in root.findall('.//wp:docPr',NS))+1));dp.set('name','Figure 10');dp.set('descr','Reconstruction timing on DNA Fountain, Organick File 16 and Song-format HT176, with exact recovery for Organick and HT176 and RS-valid counts for Fountain.')
 picture.find('.//pic:cNvPr',NS).set('name',image.name)
 rels=E.fromstring(parts['word/_rels/document.xml.rels']);rid='rId'+str(max(int(x.get('Id')[3:]) for x in rels if x.get('Id').startswith('rId'))+1);E.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',Id=rid,Type=NS['r']+'/image',Target='media/'+image.name);picture.find('.//a:blip',NS).set('{'+NS['r']+'}embed',rid)
 legend=h.paragraph('Fig. 10 | Reconstruction speed on additional DNA storage datasets. a, Elapsed reconstruction time after shared preprocessing; bars show means and error bars sample s.d. of three serial repeats, on a logarithmic time axis. b, Exact encoded-region recovery over all 607,150 Organick File 16 designs and 210,000 Song HT176 designs. c, RS-zero-syndrome droplets with unique seeds for DNA Fountain; the dashed line shows the published synthesis count, rather than independently verified sequence recovery. Graph modes share k within a dataset (21, 29 and 31, respectively), use eight counting workers and single-threaded search, and disable both Bloom filters. The DBGPS rule is emulated within the same adapters. The conventional baseline uses exact-prefix clustering, one MUSCLE alignment of up to 15 reads per cluster and weighted voting with eight workers. Conflicting valid outputs are rejected. All methods start from the same retained molecules; HT176 includes only successfully merged fragments. Table 6 reports times, output quality and speed ratios. Shared preprocessing, tuning, scoring and outer file decoding are excluded.',find('Fig. 9 |'))
 s=txt(legend);split=s.index(' a,');replace(legend,s[:split]);first=legend.find('w:r',NS);rp=E.Element(W+'rPr');rp.extend([h.element('b'),h.element('sz',val=22)]);first.insert(0,rp);h.append_text(legend,s[split:]);rp=E.Element(W+'rPr');rp.append(h.element('sz',val=22));legend.findall('w:r',NS)[1].insert(0,rp)
 index=len(body)-1 if body[-1].tag==W+'sectPr' else len(body);body.insert(index,picture);body.insert(index+1,legend)
 updates={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),'word/_rels/document.xml.rels':E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)}
 with ZipFile(DEST,'w') as z:
  for info in infos:z.writestr(info,updates.get(info.filename,parts[info.filename]))
  z.writestr('word/media/'+image.name,image.read_bytes())
 with ZipFile(DEST) as z:
  assert z.testzip() is None
  for name,data in parts.items():
   if name.startswith('word/media/'):assert z.read(name)==data
 oldroot=E.fromstring(parts['word/document.xml']);m={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
 assert [E.tostring(x) for x in oldroot.findall('.//m:oMath',m)]==[E.tostring(x) for x in root.findall('.//m:oMath',m)]
 record=dict(source=str(SOURCE),output=str(DEST),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),figures=len(root.findall('.//w:drawing',NS)),tables=len(root.findall('.//w:tbl',NS)),preserved_existing_media=True,preserved_equations=True)
 (R/'document_validation.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
if __name__=='__main__':main()
