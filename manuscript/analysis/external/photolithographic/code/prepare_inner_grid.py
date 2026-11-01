#!/usr/bin/env python3
"""Certify a safe per-configuration upper bound before expensive reassembly."""
import csv
import json
import pathlib

ROOT=pathlib.Path(__file__).resolve().parents[2]
RESULTS=ROOT/'input/antkowiak2020/results'
OUT=RESULTS/'inner_rs';OUT.mkdir(exist_ok=True)
rows=json.loads((RESULTS/'file3_graph_metrics.json').read_text())
lower_bound=max(r['correct'] for r in rows)
assert lower_bound==137662 and len(rows)==179, 'Original experiment grid changed'
certificate=[]
for r in rows:
    path=RESULTS/f"file3.k{r['k']}.c{r['cutoff']}.r{r['ratio']:.2f}.status.tsv"
    total=0;nonoverflow=0;seen=0
    with path.open() as handle:
        for s in csv.DictReader(handle,delimiter='\t'):
            assert int(s['id'])==seen;seen+=1
            if s['status']!='overflow':total+=int(s['complete_paths']);nonoverflow+=1
    assert seen==196596
    certificate.append(dict(k=r['k'],cutoff=r['cutoff'],ratio=r['ratio'],strict_correct=r['correct'],
                            complete_paths=total,nonoverflow_indices=nonoverflow,
                            dominated_by_strict_baseline=total<lower_bound))
(OUT/'file3_candidate_path_bounds.json').write_text(json.dumps(certificate,indent=2))
print('Original configurations:',len(certificate),'must rerun:',sum(not r['dominated_by_strict_baseline'] for r in certificate))
