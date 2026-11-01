#!/usr/bin/env python3
"""Adaptive residual-seed refinement; no reference identity selects queries."""
from collections import defaultdict
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2];OPT=ROOT/'optimization'
runs=json.loads((OPT/'results/inner_full_summary.json').read_text())
votes={}
for i,run in enumerate(runs):
    with open(run['output']) as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            body=row['payload'];bit=1<<(2*i+(row['direction']=='R'))
            votes[body]=votes.get(body,0)|bit
groups=defaultdict(list)
for body,mask in votes.items():groups[body[:13]].append((body,bin(mask).count('1')))
stable={key for key,values in groups.items() if len(values)==1 and values[0][1]>=3}
residual=OPT/'tmp/residual_forward.seeds.tsv';counts={'prior_stable_headers':len(stable),'included':0,'excluded':0}
with (OPT/'tmp/all.seeds.tsv').open() as source,residual.open('w') as dest:
    for line in source:
        direction,seed,support=line.split()
        if direction!='F':continue
        if seed[-14:-1] in stable:counts['excluded']+=1
        else:dest.write(line);counts['included']+=1
(OPT/'results/refinement_seed_provenance.json').write_text(json.dumps({**counts,
    'policy':'forward observed seeds whose 13 nt encoded header lacks a unanimous prediction supported by >=3 previous parameter/direction sources',
    'reference_sequences_used_to_select_queries':False},indent=2))
print(counts,flush=True)
# Load the reference scorer only after writing prediction-only query selection.
from optimize import Graph
graph=Graph(31,'inner_refine')
results=[]
for cutoff,ratios,skips in [(0,[1.4,1.6,1.8,2.1,2.25,2.55,3.0],[0,4]),(1,[1.4,1.6,1.8,2.1,2.55],[4])]:
    for ratio in ratios:
        for skip in skips:results.append(graph.run(cutoff,ratio,skip,residual))
graph.close()
(OPT/'results/refined_summary.json').write_text(json.dumps(results,indent=2))
print('Residual refinement complete',flush=True)
