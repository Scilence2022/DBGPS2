"""Audit and package the completed exploration without changing measurements."""
from pathlib import Path
import csv
import datetime
import hashlib
import json
import shutil
import zipfile

R = Path(__file__).resolve().parent
ARCHIVE = R.parent / 'DBGPS2_reconstruction_boundary_exploration.zip'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    v = json.loads((R / 'validation_summary.json').read_text())
    assert v['synthetic_graphs'] == 1368 and v['synthetic_operating_points'] == 43776
    assert v['real_operating_points'] == 16
    for field in ['all_native_sample_comparisons_passed', 'exact_bound_inequalities_passed',
                  'loss_decomposition_passed', 'admissible_search_losses_equal_overflow_plus_ambiguity']:
        assert v[field], field
    assert json.loads((R / 'input_revalidation.json').read_text())['all_match']
    for name, sha in json.loads((R / 'execution_manifest.json').read_text())['files'].items():
        assert digest(R / name) == sha, name
    for lock in ['geometric_source_lock.json', 'regular_geometry_source_lock.json']:
        for name, sha in json.loads((R / lock).read_text()).items():
            assert digest(R / name) == sha, name
    source = json.loads((R / 'source_provenance.json').read_text())
    for name, sha in source['core'].items():
        assert digest(R / 'src' / name) == sha, name
    assert digest(R / 'src' / 'original_profile.c') == source['original_profiler_sha256']

    records = [p for p in (R / 'runs').glob('*.json') if not p.name.endswith('.oracle.json')]
    assert len(records) == 1368
    for p in records:
        data = json.loads(p.read_text())
        assert len(data['grid']) == 32, p.name
        assert len(data['input_sha256']) == 64
        reference = R / 'references' / Path(data['command'][data['command'].index('-x') + 1]).name
        assert reference.is_file(), p.name
        if 'reference_sha256' in data:
            assert digest(reference) == data['reference_sha256'], p.name
    with (R / 'all_simulation_results.csv').open() as f:
        assert sum(1 for _ in csv.DictReader(f)) == 43776

    # Visual inspection was performed on all three delivered PNGs. Hash the
    # inspected versions so any later regeneration invalidates this QA record.
    figures = []
    for name in ['Fig_boundary_phase_maps', 'Fig_boundary_ratio_window', 'Fig_boundary_fragment_geometry']:
        files = {ext: dict(bytes=(R / f'{name}.{ext}').stat().st_size,
                           sha256=digest(R / f'{name}.{ext}')) for ext in ['png', 'pdf', 'svg']}
        figures.append(dict(name=name, png_visually_inspected=True, formats=files,
                            unresolved_layout_issues=[]))
    (R / 'figure_validation.json').write_text(json.dumps(dict(
        checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_passed=True, method='Visual inspection of the delivered PNGs; vector files present.',
        figures=figures), indent=2))

    reference_dir = R / 'reference_evaluation_only'
    reference_dir.mkdir(exist_ok=True)
    shutil.copy2('/Users/song/dbgps2_ht/ref.fa', reference_dir / 'HT176_ref.fa')
    validation_dir = R / 'validation_reference'
    validation_dir.mkdir(exist_ok=True)
    shutil.copy2('/Users/song/dbgps2_ht/results/full/HT176.grid.tsv', validation_dir / 'HT176.grid.tsv')
    prior = list(csv.DictReader((validation_dir / 'HT176.grid.tsv').open(), delimiter='\t'))
    real = list(csv.DictReader((R / 'real_HT176_raw.grid.tsv').open(), delimiter='\t'))
    matches = 0
    for row in real:
        old = next((p for p in prior if float(p['r']) == float(row['r']) and p['c'] == row['c']), None)
        if old is not None:
            for field in ['correct', 'wrong', 'no_path', 'overflow', 'no_valid', 'ambiguous']:
                assert row[field] == old[field], (row['r'], row['c'], field)
            matches += 1
    assert matches == 12

    allowed_dirs = {'runs', 'src', 'bin', 'references', 'reference_evaluation_only', 'validation_reference'}
    files = sorted(p for p in R.rglob('*') if p.is_file() and
                   (len(p.relative_to(R).parts) == 1 or p.relative_to(R).parts[0] in allowed_dirs) and
                   p.name != 'archive_manifest.json')
    manifest = dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    measurements=v, files={str(p.relative_to(R)): dict(bytes=p.stat().st_size,
                                                                     sha256=digest(p)) for p in files})
    manifest_path = R / 'archive_manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2))
    files.append(manifest_path)
    temp = ARCHIVE.with_suffix('.zip.tmp')
    with zipfile.ZipFile(temp, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6,
                         allowZip64=True) as z:
        for p in files:
            z.write(p, 'reconstruction_boundary/' + str(p.relative_to(R)))
    print('ARCHIVE_WRITTEN', len(files), 'files', flush=True)
    # Read every archived entry: ZIP CRC checks and independent SHA256 checks.
    with zipfile.ZipFile(temp) as z:
        assert len(z.infolist()) == len(files)
        for info in z.infolist():
            name = str(Path(info.filename).relative_to('reconstruction_boundary'))
            expected = manifest['files'].get(name)
            if name == 'archive_manifest.json':
                expected = dict(bytes=manifest_path.stat().st_size, sha256=digest(manifest_path))
            assert expected is not None
            h = hashlib.sha256()
            with z.open(info) as f:
                for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
                    h.update(block)
            assert info.file_size == expected['bytes'] and h.hexdigest() == expected['sha256'], name
    temp.replace(ARCHIVE)
    result = dict(archive=str(ARCHIVE), bytes=ARCHIVE.stat().st_size, sha256=digest(ARCHIVE),
                  files=len(files), all_entry_crc_and_sha256_checks_passed=True,
                  previous_raw_profile_points_revalidated=matches,
                  uncompressed_bytes=sum(p.stat().st_size for p in files),
                  excluded=['raw FASTQ', 'scratch', 'smoke development inputs',
                            'excluded_initial_fixed_lookahead_oracle', '__pycache__'])
    ARCHIVE.with_suffix('.validation.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
