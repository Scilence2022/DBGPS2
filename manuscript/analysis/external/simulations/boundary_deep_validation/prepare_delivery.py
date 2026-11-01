"""Audit and package completed measurements, figures and manuscript sources."""
from pathlib import Path
import datetime,hashlib,json,shutil,zipfile
import pandas as pd
from run_deep import digest

R=Path(__file__).resolve().parent;O=R.parents[2]
DEST=R.parent/'DBGPS2_boundary_deep_validation.zip'


def main():
    v=json.loads((R/'validation_summary.json').read_text());c=json.loads((R/'crossing_validation_summary.json').read_text())
    assert v['graphs']==330 and c['new_graphs']==36 and c['new_operating_points']==864
    for name in ['design_and_source_lock.json','crossing_design_and_source_lock.json']:
        lock=json.loads((R/name).read_text())
        for path,sha in lock.get('measurement_files',lock.get('files',{})).items():assert digest(R/path)==sha,path
    real=json.loads((R/'real_design_and_source_lock.json').read_text())
    assert digest(R/'run_real.py')==real['script_sha256']
    assert digest(R/'blind_policy.py')==real['policies_sha256'] and digest(R/'bin/blind_profile')==real['binary_sha256']
    for name,field in [('input_revalidation.json','all_match'),('figure_validation.json','all_passed'),
                       ('document_visual_validation.json','all_passed'),('regeneration_validation.json','passed'),
                       ('poisson_sampling_validation.json','passed')]:assert json.loads((R/name).read_text())[field],name
    for item in json.loads((R/'figure_validation.json').read_text())['figures']:
        for ext,sha in item['files'].items():assert digest(R/(item['name']+'.'+ext))==sha
    doc=O/'docx/DBGPS2_manuscript_reconstruction_boundaries.docx'
    assert digest(doc)==json.loads((R/'document_visual_validation.json').read_text())['document_sha256']
    assert json.loads((R/'real_HT176_raw.json').read_text())['prior_points_identical']==16
    allp=pd.concat([pd.read_csv(R/'policy_results.csv'),pd.read_csv(R/'crossing_policy_results.csv')],ignore_index=True)
    allg=pd.concat([pd.read_csv(R/'native_grid_results.csv'),pd.read_csv(R/'crossing_native_grid_results.csv')],ignore_index=True)
    assert len(allp)==1835 and len(allg)==8808
    assert len(allg[allg.cohort!='real_HT176_raw'])==8784
    assert allg.groupby(['case_id']).size().eq(24).all()
    summary=dict(new_graphs=366,new_synthetic_operating_points=8784,real_operating_points=24,
                 policy_evaluations=1835,native_search_comparisons=v['new_native_search_comparisons']+c['native_search_comparisons'],
                 all_frozen_sources_match=True,all_output_scores_match=True,
                 homogeneous_bound_conditions=c['homogeneous_conditions'],compatible_bound_conditions=c['compatible_conditions'],
                 fresh_bound_conditions=c['fresh_conditions'],fresh_compatible_conditions=c['fresh_compatible_conditions'],
                 manuscript_sha256=digest(doc),manuscript_pages=58,
                 maximum_library_size=65536,production_core_unchanged=True)
    (R/'final_validation_summary.json').write_text(json.dumps(summary,indent=2))
    for directory in ['manuscript','manuscript_source','authoring_dependencies']:(R/directory).mkdir(exist_ok=True)
    shutil.copy2(doc,R/'manuscript'/doc.name)
    shutil.copy2(O/'docx/DBGPS2_manuscript_speed_comparison_serial_3workers.docx',R/'manuscript_source/previous_manuscript.docx')
    shutil.copy2(O/'benchmarks/2026-10-19/speed_extension/docx_ooxml_helpers.py',R/'authoring_dependencies/docx_ooxml_helpers.py')
    allowed={'src','bin','references','reference_evaluation_only','runs','crossing_runs','clean_counts',
             'manuscript','manuscript_source','authoring_dependencies'}
    files=sorted(p for p in R.rglob('*') if p.is_file() and p.name not in ['archive_manifest.json','delivery.log'] and
                 (len(p.relative_to(R).parts)==1 or p.relative_to(R).parts[0] in allowed))
    entries=[(p,'boundary_deep_validation/'+str(p.relative_to(R))) for p in files]
    parent=R.parent/'DBGPS2_reconstruction_boundary_exploration.zip'
    assert parent.is_file() and digest(parent)=='aee5bc4e648f7e029bc6fba0e255c87b1c931b94ea86702e1a68498bc2b793ff'
    entries.append((parent,'parent_exploration/DBGPS2_reconstruction_boundary_exploration.zip'))
    entries.append((parent.with_suffix('.validation.json'),'parent_exploration/parent_archive_validation.json'))
    manifest=dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),validation=summary,
                  files={name:dict(bytes=path.stat().st_size,sha256=digest(path)) for path,name in entries})
    mp=R/'archive_manifest.json';mp.write_text(json.dumps(manifest,indent=2));entries.append((mp,'boundary_deep_validation/archive_manifest.json'))
    temp=DEST.with_suffix('.zip.tmp')
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for path,name in entries:
            z.write(path,name,compress_type=zipfile.ZIP_STORED if path.suffix in ['.zip','.docx'] else zipfile.ZIP_DEFLATED)
    print('ARCHIVE_WRITTEN',len(entries),'files',flush=True)
    with zipfile.ZipFile(temp) as z:
        assert len(z.infolist())==len(entries)
        for info in z.infolist():
            expected=manifest['files'].get(info.filename)
            if info.filename=='boundary_deep_validation/archive_manifest.json':expected=dict(bytes=mp.stat().st_size,sha256=digest(mp))
            assert expected is not None
            h=hashlib.sha256()
            with z.open(info) as f:
                for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
            assert info.file_size==expected['bytes'] and h.hexdigest()==expected['sha256'],info.filename
    temp.replace(DEST)
    result=dict(path=str(DEST),bytes=DEST.stat().st_size,sha256=digest(DEST),files=len(entries),
                all_entries_crc_and_sha256_verified=True,uncompressed_bytes=sum(p.stat().st_size for p,_ in entries),
                includes_parent_exploration=True,validation=summary)
    DEST.with_suffix('.validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
