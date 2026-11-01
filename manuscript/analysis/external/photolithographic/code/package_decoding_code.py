#!/usr/bin/env python3
"""Package per-file decoder sources with an offline build and explicit replay stages."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'input/antkowiak2020'
OUT = ROOT / 'output/antkowiak2020_reports'
STAGING = ROOT / 'tmp/report_packages'
BENCH = ROOT / 'benchmarks/antkowiak2020'

BOOTSTRAP = r'''#!/usr/bin/env python3
"""Build packaged snapshots without a Git checkout or network access."""
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'input/antkowiak2020';BENCH=ROOT/'benchmarks/antkowiak2020'
for name in ['bin','results','processed','raw']:(DATA/name).mkdir(parents=True,exist_ok=True)
def run(args):subprocess.run([str(a) for a in args],cwd=ROOT,check=True)
common=[ROOT/a for a in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']]
run(['cc','-O3','-Wall','-I.',BENCH/'graph_adapter.c',*common,'-lz','-lpthread','-o',DATA/'bin/graph_adapter'])
run(['cc','-O3',BENCH/'prepare_reads.c','-lz','-o',DATA/'bin/prepare_reads'])
for fast in [False,True]:
 run(['c++','-O3','-std=c++11',*(['-DFAST_GF14'] if fast else []),'-I'+str(BENCH),'-I'+str(DATA/'sources/portable_rs'),BENCH/'rs_driver.cpp','-o',DATA/('bin/rs_driver_fast' if fast else 'bin/rs_driver')])
run(['c++','-O3','-std=c++11',BENCH/'inner_select.cpp','-o',DATA/'bin/inner_select'])
run(['c++','-O3','-std=c++11','-I'+str(BENCH),BENCH/'test_inner_graph.cpp',BENCH/'graph_inner_rs.cpp','-o',DATA/'bin/test_inner_graph'])
run(['c++','-O3','-std=c++11',BENCH/'test_fast_gf.cpp','-I'+str(BENCH),'-I'+str(DATA/'sources/portable_rs'),'-o',DATA/'bin/test_fast_gf'])
run(['python3',BENCH/'build_inner_graph.py'])
run(['python3',BENCH/'build_block_graph.py'])
if (ROOT/'vendor/DBGPS-analyzer').exists():
 source=ROOT/'vendor/DBGPS-analyzer';dest=DATA/'bin/analyzer_tools';dest.mkdir(exist_ok=True)
 run(['cc','-O3','-Wall','-Wextra',source/'DBGPS-analyzer.c',source/'kthread.c',source/'vendor/bfc/bbf.c','-lz','-lpthread','-o',dest/'DBGPS-analyzer'])
 run(['cc','-O3','-Wall','-Wextra',source/'DBGPS-links.c',source/'kthread.c','-lz','-lpthread','-o',dest/'DBGPS-links'])
 run(['c++','-O3','-std=c++11',BENCH/'split_aligned_blocks.cpp','-lz','-o',DATA/'bin/split_aligned_blocks'])
print('PACKAGED BUILD COMPLETE')
'''

REPLAY = r'''#!/usr/bin/env python3
"""Explicit stages; no download or expensive decoding occurs without invocation."""
import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'input/antkowiak2020';BENCH=ROOT/'benchmarks/antkowiak2020'
sys.path.insert(0,str(BENCH))
PACKAGE=json.loads((ROOT/'PACKAGE.json').read_text());F=PACKAGE['file'];label='file'+str(F)
RAW={1:'I16_S2_R1_001.fastq.gz',2:'I18_S3_R1_001.fastq.gz',3:'File3.fq.gz'}[F]
def run(args,env=None):subprocess.run([str(a) for a in args],cwd=ROOT,check=True,env=env)
def py(script,*args):run([sys.executable,BENCH/script,*args])
def graph(ks,cuts=None,ratios=None):
 args=['--k',*map(str,ks)]
 if cuts:args+=['--cutoffs',cuts]
 if ratios:args+=['--ratio',ratios]
 py('run_graphs.py',label,*args,'--parallel','2')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['build','smoke','fetch','prepare','graph','msa','outer','inner','analyzer','partition','small-k','supplement']);args=ap.parse_args()
 if args.stage=='build':run([sys.executable,ROOT/'tools/bootstrap.py']);return
 if args.stage=='smoke':
  run([DATA/'bin/test_fast_gf'])
  run([DATA/'bin/test_inner_graph',ROOT/'fixtures/File1_validation.fasta',DATA/'results/fixture_inner'])
  return
 if args.stage=='fetch':
  import fetch_sources,download
  fetch_sources.main()
  download.download(next(spec for spec in download.FILES if spec[1]==RAW))
  return
 if args.stage=='prepare':
  py('prepare_references.py')
  if F in [2,3]:
   run([DATA/'bin/rs_driver_fast','decode',DATA/f'sources/File{F}_synthesized.fasta',DATA/f'results/reference_file{F}.zip',str(PACKAGE['blocks'])])
  with (DATA/f'results/{label}_reads.json').open('w') as dst:subprocess.run([str(DATA/'bin/prepare_reads'),str(DATA/'raw'/RAW),str(DATA/'processed'/label)],stdout=dst,check=True)
  return
 if args.stage=='graph':
  if F==1:
   py('run_graphs.py',label,'--parallel','2')
   graph([13,15,17,19,21,25,31],'9,19,49,99,199')
   graph([15,17,19,21,25,31],'1,3,4,5,6')
  else:graph([13,15,17,19,21,25,31],'1,3,9,19,49')
  if F==3:
   graph([15,17,19,21],'1,3,9,19','1.2,1.5,1.8,2.1,3.5')
   graph([13],'9,49,99','1.2,1.5,1.8')
   graph([18,20],'1,3,9','1.8,2.1,2.55,3,3.5,5,10')
   graph([19],'1,3,9','3,5,10')
   graph([18],'2','1.8,2.1,2.55,3.5')
  return
 if args.stage=='msa':
  py('build_muscle.py')
  ref=DATA/f'sources/File{F}_design.fasta' if F==1 else DATA/f'sources/File{F}_synthesized.fasta'
  py('msa_baseline.py',DATA/f'processed/{label}.fa.gz',ref,DATA/'bin/muscle3',DATA/f'results/{label}_msa','--workers','8')
  run([DATA/'bin/inner_select',DATA/f'results/{label}_msa.candidates.txt',str(PACKAGE['strands']),DATA/f'results/{label}_msa_inner'])
  if F==1:
   py('msa_baseline.py',DATA/'raw'/RAW,ref,DATA/'bin/muscle3',DATA/'results/file1_notebook','--workers','8','--prefix-length','14','--min-cluster','4','--gap-weight','0.5','--raw-fastq')
   run([DATA/'bin/inner_select',DATA/'results/file1_notebook.candidates.txt','16383',DATA/'results/file1_notebook_inner'])
  return
 if args.stage=='outer':
  if F==3:raise SystemExit('File3 has no verified sequencing-only full recovery. Use supplement after partition/small-k, with its explicit simulation labels.')
  name='file1.k19.c1.r2.55.strands' if F==1 else 'file2.k21.c1.r2.55.strands'
  run([DATA/'bin/rs_driver_fast','decode',DATA/'results'/name,DATA/f'results/{label}_packaged_recovered.zip',str(PACKAGE['blocks'])])
  return
 if F!=3:raise SystemExit('This advanced stage was measured only for File3.')
 if args.stage=='inner':
  py('prepare_inner_grid.py');py('run_inner_graph.py','--parallel','3')
  py('report_inner_graph.py');py('run_inner_graph.py','--remaining','--k','19','20','18','17','15','13','21','25','31','--parallel','3')
  py('report_parameter_union.py');py('validate_parameter_union.py');return
 if args.stage=='analyzer':
  import run_block12_analyzer as analyzer
  analyzer.REPO=ROOT/'vendor/DBGPS-analyzer';sys.argv=[str(BENCH/'run_block12_analyzer.py'),'--parallel','3'];analyzer.main()
  py('report_block12_analyzer.py');return
 if args.stage=='partition':
  py('split_reads_by_alignment.py');py('run_block_decode.py','--phase','controls','--parallel','2')
  py('run_block_decode.py','--phase','assigned-control','--parallel','2');py('run_block_decode.py','--phase','blocks','--parallel','3');py('report_block_decode.py');return
 if args.stage=='small-k':
  for stage in ['initial','fine','fine14']:py('run_block12_small_k.py','--stage',stage)
  py('report_block12_small_k.py');return
 if args.stage=='supplement':py('simulate_file3_supplement.py');return
if __name__=='__main__':main()
'''


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    common = ['README.md','fetch_sources.py','download.py','build.py','build_muscle.py','prepare_reads.c','prepare_references.py',
              'graph_adapter.c','run_graphs.py','msa_baseline.py','inner_select.cpp','rs_driver.cpp','operators_compat.hpp','legacy_random.hpp',
              'fast_gf14.hpp','test_fast_gf.cpp','graph_inner_rs.cpp','graph_inner_rs.h','test_inner_graph.cpp','build_inner_graph.py','build_block_graph.py','report.py']
    core = ['DBGPS2.c','dbgps_codec.c','dbgps_codec.h','dbgps_redundancy.c','dbgps_redundancy.h','dbgps_strand.c','dbgps_strand.h',
            'dbgps_bloom.h','dbgps_bloom_graph.h','kthread.c','kthread.h','ketopt.h','khashl.h','kseq.h','LICENSE']
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    manifest = []
    for number in [1,2,3]:
        package = STAGING / f'File{number}_code'
        if package.exists(): shutil.rmtree(package)
        package.mkdir(parents=True)
        files = common if number < 3 else [p.name for p in BENCH.iterdir() if p.is_file() and p.suffix in ['.py','.c','.cpp','.h','.hpp','.md'] and p.name not in ['write_full_reports.py','package_decoding_code.py']]
        def copy(source, relative):
            destination = package / relative; destination.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source,destination)
        for name in core: copy(ROOT/name,name)
        for name in files: copy(BENCH/name,'benchmarks/antkowiak2020/'+name)
        for source in (DATA/'sources/portable_rs').glob('*.hpp'): copy(source,'input/antkowiak2020/sources/portable_rs/'+source.name)
        for name in ['outer_probe.cpp','outer_spectrum.cpp','rng.cpp','test_legacy.cpp']:
            copy(DATA/'results'/name,'diagnostic_sources/'+name)
        copy(DATA/'sources/muscle3.8.31_src.tar.gz','input/antkowiak2020/sources/muscle3.8.31_src.tar.gz')
        copy(DATA/'sources/noisy_dna_data_storage/LICENSE','third_party_licenses/noisy_dna_data_storage.LICENSE')
        # The upstream RS checkout has no standalone LICENSE file; preserve its provenance and README.
        copy(DATA/'sources/dna_rs_coding/README.md','third_party_licenses/dna_rs_coding.README.md')
        fixture=package/'fixtures/File1_validation.fasta';fixture.parent.mkdir(exist_ok=True)
        sequences=[line for line in (DATA/'sources/File1_design.fasta').read_text().splitlines() if line and not line.startswith('>')]
        assert len(sequences)==16383 and all(len(s)==60 for s in sequences)
        fixture.write_text(''.join(f'>{i:06d}\n{s}\n' for i,s in enumerate(sequences)))
        tools=package/'tools';tools.mkdir();(tools/'bootstrap.py').write_text(BOOTSTRAP);(tools/'reproduce.py').write_text(REPLAY)
        (package/'requirements-analysis.txt').write_text('numpy>=2.0,<3\nmatplotlib>=3.9,<4\n')
        if number==3:
            analyzer=Path('/Users/song/.codex/worktrees/dabb/DBGPS-analyzer')
            provenance=json.loads((DATA/'results/analyzer_block12/provenance.json').read_text())
            for name, expected in provenance['source_sha256'].items():
                assert sha(analyzer/name)==expected,name
                copy(analyzer/name,'vendor/DBGPS-analyzer/'+name)
            for name in ['ketopt.h','kthread.h','LICENSE','vendor/bfc/LICENSE.txt']:
                copy(analyzer/name,'vendor/DBGPS-analyzer/'+name)
            # This is an auxiliary input read by the block report, not decoder evidence fabricated on a fresh replay.
            copy(DATA/'results/block_decode/raw_ancillary_scope.json','recorded_evidence/raw_ancillary_scope.json')
        meta=dict(file=number,strands={1:16383,2:49149,3:196596}[number],blocks={1:1,2:3,3:12}[number],
                  dbgps2_commit=commit,adapter='Experimental internal-index adapter; not stock CLI support',
                  rs_headers_revision='8e7e3b5',legacy_sampler=True,raw_reads_included=False,
                  dependencies='C compiler, C++11, zlib, pthread, Python 3.12 recommended; NumPy/Matplotlib for advanced reporting; minimap2 2.31-r1302 for File3 partitioning.',
                  packaging_changes='Added tools/bootstrap.py and tools/reproduce.py; copied all original decoder scripts unchanged. Vendored prepared portable RS headers avoid a Git/network requirement at build time. Analyzer wrapper relocates its provenance path in memory only. The validation FASTA preserves all File1 sequences with normalized headers/blank lines; the test mutates its first 100 words and uses the full domain for decision controls.')
        (package/'PACKAGE.json').write_text(json.dumps(meta,indent=2))
        readme=f'''# File {number} decoding code / 解码代码包

This is the recorded experimental adapter, not the stock DBGPS2 CLI.
本包保留原解码脚本，新增离线构建与显式分阶段入口；未加入原始 FASTQ。

Run from this package root / 从包根目录运行:

```sh
python3 tools/reproduce.py build
python3 tools/reproduce.py smoke
# Network and full-data computation occur only in explicitly requested stages:
python3 tools/reproduce.py fetch
python3 tools/reproduce.py prepare
python3 tools/reproduce.py graph
python3 tools/reproduce.py msa
'''+('python3 tools/reproduce.py outer\n' if number<3 else '''python3 -m pip install -r requirements-analysis.txt
python3 tools/reproduce.py inner
python3 tools/reproduce.py analyzer
python3 tools/reproduce.py partition
python3 tools/reproduce.py small-k
python3 tools/reproduce.py supplement
''')+'''```

`build` is offline and uses the exact prepared RS-header snapshot from the experiment.
`smoke` runs field-equivalence and 126000 single-symbol correction checks on 100 public File1 codewords.
`fetch` retrieves pinned public sources and ONLY this file's complete FASTQ, checking its MD5. Small design/source files for all three datasets are fetched by the original shared fetcher.
`prepare` reconstructs the synthesized evaluation references and processes this input's reads.
`graph` replays the exact final original grid: 79/35/179 settings for File1/2/3, respectively.
`msa` also builds the pinned MUSCLE 3.8.31 archive; it follows the original fixed-seed baseline, plus the notebook variant for File1.
`outer` uses the best actual strict output for File1 or File2. A ZIP's existence alone does not establish correctness: compare bytes and CRC to the validated baseline recorded in the report.

For File3 advanced stages, all original 179 pooled status/metric files are prerequisites to inner analysis. Partitioning requires `minimap2` on PATH, version 2.31-r1302 for the recorded configuration. The preserved mapper command uses `samtools` only if the original script requests it; check that tool on PATH as applicable.
`analyzer` uses the bundled source snapshot, relocating only the provenance path in memory.
`partition` completes the 319-setting core including the matched same-read re-pooled control; expensive optional unfiltered raw controls are not enabled automatically.
`small-k` requires completed broad block reads and reports 174 added settings plus 179 earlier settings.
`supplement` requires the block and small-k reports. It deliberately injects design strands and creates files named `simulated_recovered.zip`; it is NOT sequencing-only recovery.

For File1 there is one block, so cross-block reads splitting is inapplicable. File2 was not measured after reads splitting. File3 is the only dataset with this task's alignment-routing experiments.

The original `report.py` combines ALL THREE datasets and is retained for provenance. Use the delivered per-file reports for the frozen results. Original scripts containing absolute local example paths are retained unchanged; use the new wrappers for portable invocation. `/usr/bin/time -l` used by historical run scripts is macOS-specific; Linux users must adapt only the resource-recording wrapper, without changing algorithms.

No raw FASTQ, compiled binaries, caches, manuscript files or credentials are included. Sources are accompanied by the available original licenses and provenance. The RS upstream checkout has no standalone LICENSE file; its README/source comments are preserved without inventing a license.

See PACKAGE.json for versions and SOURCE_FILES.json for SHA-256 file provenance. `requirements-analysis.txt` declares compatible ranges for reporting, not a bit-for-bit historical Python environment lock.
'''
        (package/'README.md').write_text(readme)
        source_manifest={str(p.relative_to(package)):sha(p) for p in sorted(package.rglob('*')) if p.is_file()}
        (package/'SOURCE_FILES.json').write_text(json.dumps(source_manifest,indent=2))
        destination=OUT/'code'/f'File{number}_code.zip';destination.parent.mkdir(exist_ok=True)
        with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
            for source in sorted(package.rglob('*')):
                if source.is_file():archive.write(source,str(source.relative_to(STAGING)))
        manifest.append(dict(file=number,archive=str(destination),files=len(source_manifest)+1,bytes=destination.stat().st_size,sha256=sha(destination)))
        print('PACKAGED',number,'files',len(source_manifest)+1,'bytes',destination.stat().st_size,flush=True)
    (OUT/'code/package_manifest.json').write_text(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
