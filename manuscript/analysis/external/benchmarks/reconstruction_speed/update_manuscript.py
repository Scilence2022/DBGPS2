"""Patch the integrated manuscript without rebuilding its existing OOXML."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

from lxml import etree as E
from PIL import Image

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[2]/'docx/DBGPS2_manuscript_with_figures_Gimpel2024.docx'
DESTINATION=SOURCE.with_name('DBGPS2_manuscript_speed_comparison.docx')
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing',
    'pic':'http://schemas.openxmlformats.org/drawingml/2006/picture',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
W='{'+NS['w']+'}'


def element(tag,**attrs):
    return E.Element(W+tag,{W+k:str(v) for k,v in attrs.items()})


def text(node):
    return ''.join(node.itertext())


def replace_text(paragraph,value):
    for child in list(paragraph):
        if child.tag!=W+'pPr': paragraph.remove(child)
    run=E.SubElement(paragraph,W+'r');t=E.SubElement(run,W+'t')
    t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=value


def append_text(paragraph,value):
    run=E.SubElement(paragraph,W+'r');t=E.SubElement(run,W+'t')
    t.set('{http://www.w3.org/XML/1998/namespace}space','preserve');t.text=value


def replace_run_text(paragraph,old,new):
    matches=[t for t in paragraph.findall('.//w:t',NS) if old in (t.text or '')]
    assert len(matches)==1,(old,len(matches))
    matches[0].text=matches[0].text.replace(old,new)


def paragraph(value,template):
    result=deepcopy(template);replace_text(result,value)
    first_run=template.find('w:r',NS)
    properties=first_run.find('w:rPr',NS) if first_run is not None else None
    if properties is not None:result.find('w:r',NS).insert(0,deepcopy(properties))
    return result


def table(rows):
    widths=[900,2900,1900,1400,900,1360]
    result=element('tbl');properties=E.SubElement(result,W+'tblPr')
    properties.extend([element('tblW',w=sum(widths),type='dxa'),element('jc',val='center'),element('tblLayout',type='fixed')])
    borders=E.SubElement(properties,W+'tblBorders')
    for side in ['top','left','bottom','right','insideH','insideV']:
        borders.append(element(side,val='single',sz=4,color='D9D9D9'))
    margins=E.SubElement(properties,W+'tblCellMar')
    for side in ['top','bottom']:margins.append(element(side,w=100,type='dxa'))
    for side in ['left','right']:margins.append(element(side,w=100,type='dxa'))
    grid=E.SubElement(result,W+'tblGrid')
    for width in widths:grid.append(element('gridCol',w=width))
    for row_index,values in enumerate(rows):
        row=E.SubElement(result,W+'tr');row_properties=E.SubElement(row,W+'trPr')
        row_properties.append(element('cantSplit'))
        if row_index==0:row_properties.append(element('tblHeader'))
        for column,(value,width) in enumerate(zip(values,widths)):
            cell=E.SubElement(row,W+'tc');cell_properties=E.SubElement(cell,W+'tcPr')
            cell_properties.extend([element('tcW',w=width,type='dxa'),element('vAlign',val='center')])
            if row_index==0:cell_properties.append(element('shd',fill='DCE6F1',val='clear'))
            elif row_index%2==0:cell_properties.append(element('shd',fill='F5F7F9',val='clear'))
            p=E.SubElement(cell,W+'p');ppr=E.SubElement(p,W+'pPr')
            ppr.extend([element('spacing',before=0,after=0,line=240,lineRule='auto'),element('jc',val='left' if column==1 else 'center')])
            run=E.SubElement(p,W+'r');rpr=E.SubElement(run,W+'rPr')
            rpr.extend([element('rFonts',ascii='Times New Roman',hAnsi='Times New Roman'),element('sz',val=20),element('color',val='000000')])
            if row_index==0:rpr.append(element('b'))
            t=E.SubElement(run,W+'t');t.text=str(value)
    return result


def main():
    summary=json.loads((ROOT/'summary.json').read_text())
    get=lambda d,m:next(x for x in summary if x['dataset']==d and x['method']==m)
    a,b=get('file1','dbgps2'),get('file2','dbgps2')
    ga,gb=get('file1','dbgps_rule'),get('file2','dbgps_rule')
    ma,mb=get('file1','cluster_msa_vote'),get('file2','cluster_msa_vote')
    f=lambda row:f'{row["mean_seconds"]:.1f} ± {row["sd_seconds"]:.1f}'
    pct=lambda row:f'{100*row["exact_recovery"]:.2f}%'
    with ZipFile(SOURCE) as z:
        original={n:z.read(n) for n in z.namelist()};infos=z.infolist()
    root=E.fromstring(original['word/document.xml']);body=root.find('w:body',NS)
    ps=root.findall('w:body/w:p',NS)
    assert len(ps)==196 and text(ps[66])=='Discussion', 'Source manuscript structure changed'
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='38aa87bf80ffdb97a1d5e1d1a24f30945d523ce44b4ef75c9ff6380ab5558add', 'Source was edited since benchmark preparation'
    append_text(ps[11], ' Controlled repeated timings additionally compare both graph pruning strategies with prefix clustering, multiple alignment and weighted majority voting on photolithographic Files 1 and 2.')
    replace_run_text(ps[48], 'This evaluation tests the adapted reconstruction pipeline; it does not isolate the contribution of the ratio bound on these data.', 'The three-file evaluation tests the adapted reconstruction pipeline; a separate controlled comparison on Files 1 and 2 evaluates the two graph pruning strategies and their computational cost.')
    replace_run_text(ps[102], 'All reported timing runs', 'The thermal-pool timing runs')
    replace_run_text(ps[130], 'Only the controlled thermal-pool experiments were used for the matched DBGPS-versus-DBGPS2 resource and algorithm comparisons.', 'The controlled thermal-pool experiments isolate the ratio rule; the File 1/2 benchmark compares resource use at dataset-tuned operating points.')
    append_text(ps[132], ' Figure 9 and Table 5 use the new serial reconstruction benchmark; source tables retain all three repeats and stage timings for each method and input.')
    append_text(ps[139], ' The accompanying reconstruction_speed package contains the additional timing harness, source and executable hashes, pilot settings, exact commands, repeated timing logs, output hashes and figure-generation script.')

    heading=paragraph('Controlled speed comparison with clustering and multiple alignment',ps[58])
    p1=paragraph('We compared both graph pruning strategies with the prefix-clustering, MUSCLE and weighted-voting pipeline on the same processed File 1 and File 2 inputs, comprising 14,299,688 and 33,404,713 reads, respectively (Fig. 9 and Table 5). Each operating point was timed three times on the same workstation, with runs executed serially and method order rotated. Timings cover compressed-input reading through strict strand output, including inner-code selection. Shared preprocessing, parameter selection, reference scoring and outer decoding are excluded. The DBGPS comparator emulates its absolute-cutoff rule within the same format adapter rather than using the original published executable.',ps[65])
    p2=paragraph(f'For File 1, clustering, MSA and voting required {f(ma)} s (mean ± s.d., n = 3), compared with {f(ga)} s for the DBGPS rule and {f(a)} s for DBGPS2. The corresponding File 2 times were {f(mb)}, {f(gb)} and {f(b)} s. Relative to the conventional baseline, the DBGPS rule was {ga["speedup_vs_msa"]:.1f}- and {gb["speedup_vs_msa"]:.1f}-fold faster, and DBGPS2 was {a["speedup_vs_msa"]:.1f}- and {b["speedup_vs_msa"]:.1f}-fold faster. These ratios compare complete reconstruction pipelines at the reported operating points, rather than path-search time alone.',ps[65])
    p3=paragraph(f'DBGPS2 exactly reconstructed {pct(a)} and {pct(b)} of designs, compared with {pct(ga)} and {pct(gb)} for the DBGPS rule and {pct(ma)} and {pct(mb)} for the conventional baseline. Table 5 reports incorrect outputs alongside recovery. All three repeats of a method returned identical strict output sets. These timings compare recovery and cost at the selected operating points, rather than equal-recovery or full-file-decoding endpoints.',ps[65])
    append_text(p3, f' Relative to the best tested absolute-cutoff settings, DBGPS2 increased exact recovery by {100*(a["exact_recovery"]-ga["exact_recovery"]):.2f} and {100*(b["exact_recovery"]-gb["exact_recovery"]):.2f} percentage points, using {a["mean_seconds"]/ga["mean_seconds"]:.2f} and {b["mean_seconds"]/gb["mean_seconds"]:.2f} times the total reconstruction time, respectively.')
    stage=paragraph(f'Clustering required {ma["cluster_seconds_mean"]:.1f} and {mb["cluster_seconds_mean"]:.1f} s in Files 1 and 2; the parallel MSA/voting stage required {ma["msa_vote_seconds_mean"]:.1f} and {mb["msa_vote_seconds_mean"]:.1f} s, accounting for {100*ma["msa_vote_seconds_mean"]/ma["mean_seconds"]:.1f}% and {100*mb["msa_vote_seconds_mean"]/mb["mean_seconds"]:.1f}% of total baseline time. For DBGPS2, graph counting required {a["count_seconds_mean"]:.1f} and {b["count_seconds_mean"]:.1f} s, seed enumeration required {a["seed_seconds_mean"]:.1f} and {b["seed_seconds_mean"]:.1f} s, and path search required {a["search_seconds_mean"]:.1f} and {b["search_seconds_mean"]:.1f} s. The gain over this conventional implementation therefore concerns the complete reconstruction workflow, including its repeated alignment calls.',ps[65])
    caption=paragraph('Table 5 Reconstruction time and exact strand outputs on matched photolithographic inputs',ps[45])
    rows=[['Input','Method','Time (s)','Exact (%)','Wrong','Speedup']]
    for d in ['file1','file2']:
        for m,label in [('cluster_msa_vote','Clustering + MSA + vote'),('dbgps_rule','DBGPS rule (adapted)'),('dbgps2','DBGPS2 (adapted)')]:
            row=get(d,m)
            rows.append([d.replace('file','File '),label,f(row),f'{100*row["exact_recovery"]:.2f}',row['incorrect'],f'{row["speedup_vs_msa"]:.1f}×'])
    tbl=table(rows)
    footnote=paragraph('Time is mean ± sample s.d. of three serial computational repeats. Exact recovery divides the number of exact 60-nt outputs by all designs; wrong counts include every erroneous strict output. Speedup divides the mean baseline time by the mean method time. File 1 uses k = 19 and File 2 k = 21 for both graph modes. DBGPS cutoffs 5 and 3 were selected from an eight-cutoff pilot grid; DBGPS2 uses c = 1 and r = 2.55. Bloom filters are disabled. MSA candidate coverage before unique selection is reported separately in the preceding photolithographic analysis.',ps[46])
    index=body.index(ps[66])
    for node in [heading,p1,p2,stage,p3,caption,tbl,footnote]:
        body.insert(index,node);index+=1

    discussion=paragraph('The controlled reconstruction benchmark extends the resource analysis beyond comparisons between graph rules. Both graph modes avoid running a separate multiple alignment for each retained prefix cluster, and are faster than the tested clustering/MSA/voting pipeline while recovering more exact strands. This implementation-level benefit should be distinguished from the incremental contribution of the ratio constraint. Because the absolute cutoff was selected separately for each pruning strategy, their end-to-end timing difference includes a recovery–search trade-off; it need not match the speed gain measured at a shared cutoff in the thermal pool. The new experiments cover two photolithographic inputs and do not establish equivalent speedups for fragmented thermal libraries, alternative clustering implementations, or complete outer decoding.',ps[73])
    append_text(discussion, ' The baseline launches MUSCLE separately for each cluster or sampled subset, so its measured cost includes process startup and communication; an in-process alignment implementation could have a different run time.')
    body.insert(body.index(ps[73]),discussion)

    mh=paragraph('Controlled reconstruction timing against clustering and multiple alignment',ps[101])
    mp1=paragraph('For the new comparison, all methods consumed the same existing gzipped FASTA inputs for File 1 and File 2, without reference-assisted partitioning. The graph adapter used eight-threaded exact counting, observed internal-index seeds, single-threaded left/right extension, the 1,050-live-path limit and zero-syndrome inner-RS validation. Both Bloom filters were disabled. The two graph modes used the same k within an input (19 for File 1; 21 for File 2). DBGPS was emulated with r = 1,000,000. A pilot scan of c ∈ {1, 3, 5, 8, 10, 15, 20, 30} selected the cutoff with the highest exact recovery for this mode, giving c = 5 and 3. DBGPS2 retained the previously used c = 1 and r = 2.55. Pilot selection was retrospective and excluded from the repeated timings; these are fixed dataset-tuned operating points, not an independent test of parameter selection.',ps[103])
    mp2=paragraph('The prefix baseline used these parameters (16-nt key, minimum cluster size five, up to 15 reads per alignment, five fixed-seed subsets for larger clusters, MUSCLE 3.8.31 and gap-vote weight 0.4). Reconstruction did not access the known designs. Candidate generation was followed by the existing inner-RS selector; only its strict zero-syndrome output was scored. Eight MSA workers were used, with numerical-library threading capped at one. The selector also generated its corrected-output file, whose computation remains included in elapsed time. Input decompression, clustering, subprocess communication, alignment, voting, validation and output writing were included. The MSA/voting stage was measured jointly as parallel wall time; worker elapsed-time sums were retained only as diagnostics.',ps[103])
    mp3=paragraph('Three repeats per method and input were run sequentially with order rotated within datasets. Elapsed time was measured around /usr/bin/time -l; the timer ended after the reconstruction process exited and before exact-reference scoring. Background user applications remained present, no cache was forcibly dropped, and pre-run process snapshots were retained. The counter uses eight hash-insertion workers plus pipelined reading and k-mer extraction; the baseline uses eight alignment workers plus coordinator overhead. Configured worker counts therefore match, but total active thread counts and CPU utilisation need not. These are repeated workstation timings rather than dedicated-node measurements with controlled cold caches or CPU affinity. Mean times, sample standard deviations, stage timings and SHA256 hashes of strict outputs were recorded. Downloads, common preprocessing, parameter selection and outer file decoding were excluded consistently. The eight-worker photolithographic comparison is separate from the three-threaded thermal-pool resource measurements in Fig. 7.',ps[103])
    index=body.index(ps[104])
    for node in [mh,mp1,mp2,mp3]:body.insert(index,node);index+=1

    # Preserve the first eight figures and append a native drawing with a new relation.
    picture=deepcopy(ps[194]);image=ROOT/'Fig9_reconstruction_speed.png'
    with Image.open(image) as im:width,height=im.size
    cx=5943600;cy=round(cx*height/width)
    for ext in picture.findall('.//wp:extent',NS)+picture.findall('.//a:xfrm/a:ext',NS):
        ext.set('cx',str(cx));ext.set('cy',str(cy))
    docpr=picture.find('.//wp:docPr',NS)
    maximum=max(int(p.get('id')) for p in root.findall('.//wp:docPr',NS))
    docpr.set('id',str(maximum+1));docpr.set('name','Figure 9')
    docpr.set('descr','Two panels compare repeated total reconstruction times and exact strand recovery for a prefix-clustering MUSCLE voting pipeline, the adapted DBGPS rule and adapted DBGPS2 on photolithographic Files 1 and 2.')
    picture.find('.//pic:cNvPr',NS).set('name',image.name)
    rels=E.fromstring(original['word/_rels/document.xml.rels'])
    relns='{http://schemas.openxmlformats.org/package/2006/relationships}'
    rid='rId'+str(max(int(x.get('Id')[3:]) for x in rels if x.get('Id').startswith('rId'))+1)
    E.SubElement(rels,relns+'Relationship',Id=rid,Type=NS['r']+'/image',Target='media/Fig9_reconstruction_speed.png')
    picture.find('.//a:blip',NS).set('{'+NS['r']+'}embed',rid)
    legend=paragraph('Fig. 9 | Controlled reconstruction timing and recovery on two photolithographic libraries. a, Total reconstruction time from the same processed gzipped FASTA input to strict strand outputs, including clustering/MSA/voting or graph counting/seed enumeration/search, inner-code decisions and output writing. Bars show mean elapsed time and error bars show sample s.d. from three serial repeats with rotated method order; the time axis is logarithmic. b, Exact 60-nt strand recovery divided by all designed indices; outputs were identical across repeats. Both graph modes use k = 19 for File 1 and k = 21 for File 2, with eight-threaded counting and single-threaded search. The DBGPS rule is emulated in the same adapter with r = 1,000,000 and pilot-selected cutoffs 5 and 3; DBGPS2 uses c = 1 and r = 2.55. Both Bloom filters are disabled. The baseline uses exact prefix16 clustering, MUSCLE 3.8.31 and weighted majority voting with eight alignment workers. Table 5 reports incorrect outputs and mean-time ratios. Common preprocessing, tuning, reference scoring and outer file decoding are excluded. These are implementation-specific operating-point comparisons, not equal-recovery benchmarks or measurements on the thermal-aging samples.',ps[195])
    legend_text=text(legend);split=legend_text.index(' a,')
    replace_text(legend,legend_text[:split])
    first=legend.find('w:r',NS);first_properties=E.Element(W+'rPr')
    first_properties.extend([element('b'),element('sz',val=22)])
    first.insert(0,first_properties)
    append_text(legend,legend_text[split:])
    second=legend.findall('w:r',NS)[1];second_properties=E.Element(W+'rPr')
    second_properties.append(element('sz',val=22));second.insert(0,second_properties)
    index=len(body)-1 if body[-1].tag==W+'sectPr' else len(body)
    body.insert(index,picture);body.insert(index+1,legend)
    updated={'word/document.xml':E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True),
             'word/_rels/document.xml.rels':E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)}
    with ZipFile(DESTINATION,'w') as z:
        for info in infos:z.writestr(info,updated.get(info.filename,original[info.filename]))
        z.writestr('word/media/Fig9_reconstruction_speed.png',image.read_bytes())
    with ZipFile(DESTINATION) as z:assert z.testzip() is None
    report={'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'output_sha256':hashlib.sha256(DESTINATION.read_bytes()).hexdigest(),
            'output':str(DESTINATION),'figures':len(root.findall('.//w:drawing',NS)),
            'tables':len(root.findall('.//w:tbl',NS)),
            'paragraphs':len(root.findall('w:body/w:p',NS))}
    (ROOT/'document_validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
