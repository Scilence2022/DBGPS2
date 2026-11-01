"""Serial, rotated-order timing with post-run exact-sequence scoring."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
DATA = Path('/Users/song/.codex/worktrees/9a75/DBGPS2/input/antkowiak2020')
DATASETS = {'file1': dict(n=16383, k=19, ref='File1_design.fasta'),
            'file2': dict(n=49149, k=21, ref='File2_synthesized.fasta')}


def truth(label):
    return [s.strip() for s in (DATA/'sources'/DATASETS[label]['ref']).read_text().splitlines()
            if s.strip() and not s.startswith('>')]


def score(path, target):
    values = dict((int(line.split('\t')[0]), line.split('\t')[1])
                  for line in path.read_text().splitlines() if line)
    correct = sum(0 <= i < len(target) and s == target[i] for i, s in values.items())
    return dict(emitted=len(values), correct=correct, incorrect=len(values)-correct,
                exact_recovery=correct/len(target),
                output_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def pilot_selection():
    chosen = {}
    for label, config in DATASETS.items():
        target = truth(label)
        rows = []
        for path in (ROOT/'pilot').glob(label+'.k*.strands'):
            cutoff = int(re.search(r'\.c(\d+)\.', path.name)[1])
            rows.append(dict(cutoff=cutoff, **score(path, target)))
        assert rows, label
        best = max(rows, key=lambda x: (x['correct'], -x['incorrect'], -x['cutoff']))
        chosen[label] = dict(dbgps_cutoff=best['cutoff'], dbgps2_cutoff=1,
                             pilot_scores=sorted(rows, key=lambda x:x['cutoff']))
    (ROOT/'selected_settings.json').write_text(json.dumps(chosen, indent=2))
    return chosen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--select-only', action='store_true')
    args = parser.parse_args()
    chosen = pilot_selection()
    if args.select_only:
        print(json.dumps(chosen, indent=2)); return
    env = os.environ.copy()
    env.update(ANTK_LAYOUT='internal', ANTK_BLOOM='0', OMP_NUM_THREADS='1',
               OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    metadata = dict(started=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                    platform=platform.platform(), python=sys.version,
                    cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),
                    cores=16, ram_bytes=137438953472, max_workers=8,
                    input_format='shared preprocessed gzipped FASTA; preprocessing excluded for every method',
                    graph_baseline='DBGPS absolute-cutoff rule emulated in the same format adapter; not the 2022 binary',
                    bloom_filters=False, repeats=3,
                    order='Latin rotation within dataset, all measured processes serial',
                    tuning='parameters fixed before repeated timing; DBGPS cutoff selected retrospectively by pilot recovery',
                    references='accessed only by post-run evaluation, never reconstruction')
    (ROOT/'benchmark_metadata.json').write_text(json.dumps(metadata, indent=2))
    records=[]
    for repeat in range(1,4):
        for dataset_index,(label,config) in enumerate(DATASETS.items()):
            methods=['dbgps_rule','dbgps2','cluster_msa_vote']
            shift=(repeat-1+dataset_index)%3
            methods=methods[shift:]+methods[:shift]
            for order,method in enumerate(methods):
                name=f'{label}_{method}_rep{repeat}'
                prefix=ROOT/'runs'/name
                result_path=Path(str(prefix)+'.result.json')
                if result_path.exists():
                    records.append(json.loads(result_path.read_text()));continue
                reads=DATA/'processed'/f'{label}.fa.gz'
                if method=='cluster_msa_vote':
                    command=[sys.executable,str(ROOT/'msa_pipeline.py'),str(reads),str(config['n']),str(prefix),'--workers','8']
                    strands=Path(str(prefix)+'.strict.strands')
                else:
                    ratio='1000000' if method=='dbgps_rule' else '2.55'
                    cutoff=chosen[label]['dbgps_cutoff' if method=='dbgps_rule' else 'dbgps2_cutoff']
                    command=[str(ROOT/'bin/graph_adapter'),str(reads),str(config['n']),str(config['k']),str(prefix),'-',ratio,str(cutoff)]
                    strands=Path(str(prefix)+f'.k{config["k"]}.c{cutoff}.r{float(ratio):.2f}.strands')
                cpu_before=subprocess.check_output(['ps','-Ao','pid,pcpu,comm','-r'],text=True)
                Path(str(prefix)+'.load_before.txt').write_text(cpu_before)
                print('START',name,flush=True)
                start=time.perf_counter()
                with open(str(prefix)+'.stdout','w') as stdout,open(str(prefix)+'.stderr','w') as stderr:
                    completed=subprocess.run(['/usr/bin/time','-l']+command,stdout=stdout,stderr=stderr,env=env)
                seconds=time.perf_counter()-start
                assert completed.returncode==0,(name,completed.returncode)
                if method=='cluster_msa_vote':
                    stages=json.loads(Path(str(prefix)+'.stages.json').read_text())
                    stage_values=dict(cluster_seconds=stages['cluster_seconds'],msa_vote_seconds=stages['msa_vote_seconds'],validation_seconds=stages['validation_seconds'])
                    stage_values['other_seconds']=max(0,seconds-stages['pipeline_seconds'])
                else:
                    log=Path(str(prefix)+'.stderr').read_text()
                    count=float(re.search(r'COUNT k=\d+ wall=([\d.]+)',log)[1])
                    seed=float(re.search(r'wall_since_start=([\d.]+)',log)[1])-count
                    search=float(re.search(r'RESULT .* wall=([\d.]+)',log)[1])
                    stage_values=dict(count_seconds=count,seed_seconds=seed,search_seconds=search,
                                      other_seconds=max(0,seconds-count-seed-search))
                # Truth is loaded after the subprocess and timer have finished.
                record=dict(dataset=label,method=method,repeat=repeat,order=order,
                            wall_seconds=seconds,command=command,exit_code=completed.returncode,
                            k=config['k'] if method!='cluster_msa_vote' else None,
                            cutoff=cutoff if method!='cluster_msa_vote' else None,
                            ratio=float(ratio) if method!='cluster_msa_vote' else None,
                            stages=stage_values,**score(strands,truth(label)))
                result_path.write_text(json.dumps(record,indent=2));records.append(record)
                (ROOT/'results.json').write_text(json.dumps(records,indent=2))
                print('DONE',name,'seconds',round(seconds,3),'correct',record['correct'],'incorrect',record['incorrect'],flush=True)
    print('ALL_DONE',flush=True)


if __name__=='__main__':
    main()
