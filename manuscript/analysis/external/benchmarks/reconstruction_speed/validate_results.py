"""Check repeat completeness, outputs, historical compatibility and source pins."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
OLD=Path('/Users/song/.codex/worktrees/9a75/DBGPS2/input/antkowiak2020/results')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    records=json.loads((ROOT/'results.json').read_text())
    assert len(records)==18
    for dataset in ['file1','file2']:
        for method in ['dbgps_rule','dbgps2','cluster_msa_vote']:
            rows=[r for r in records if r['dataset']==dataset and r['method']==method]
            assert sorted(r['repeat'] for r in rows)==[1,2,3]
            assert len({r['output_sha256'] for r in rows})==1
            assert all(r['wall_seconds']>0 and r['exit_code']==0 for r in rows)
            assert all(abs(sum(r['stages'].values())-r['wall_seconds'])<.01 for r in rows)
    compatibility={}
    for dataset in ['file1','file2']:
        baseline=ROOT/'runs'/f'{dataset}_cluster_msa_vote_rep1'
        stages=json.loads(Path(str(baseline)+'.stages.json').read_text())
        historical=json.loads((OLD/f'{dataset}_msa.json').read_text())
        for key in ['reads','clusters','candidates_60nt']:
            current_key='eligible_clusters' if key=='clusters' else key
            assert stages[current_key]==historical[key],(dataset,key)
        strict=Path(str(baseline)+'.strict.strands')
        old=OLD/f'{dataset}_msa_inner.strict.strands'
        assert strict.read_bytes()==old.read_bytes(),dataset
        compatibility[dataset]=dict(strict_output_identical_to_historical=True,
                                     strict_sha256=digest(strict),
                                     reads=stages['reads'],clusters=stages['eligible_clusters'],
                                     candidates_60nt=stages['candidates_60nt'])
        k=19 if dataset=='file1' else 21
        graph=ROOT/'runs'/f'{dataset}_dbgps2_rep1.k{k}.c1.r2.55.strands'
        historical_graph=OLD/f'{dataset}.k{k}.c1.r2.55.strands'
        assert graph.read_bytes()==historical_graph.read_bytes(),(dataset,'graph historical compatibility')
        compatibility[dataset]['dbgps2_strict_output_identical_to_historical']=True
    pins=json.loads((ROOT/'source_sha256.json').read_text())
    assert all(digest(ROOT/p)==h for p,h in pins.items())
    report=dict(measured_runs=18,all_exited_successfully=True,
                identical_strict_outputs_across_repeats=True,
                historical_baseline_compatibility=compatibility,
                source_and_executable_pins_valid=True)
    (ROOT/'results_validation.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
