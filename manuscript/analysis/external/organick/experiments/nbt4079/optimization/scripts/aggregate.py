#!/usr/bin/env python3
"""Reference-free multi-parameter voting; references select global rules only."""
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import time
import numpy as np
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('--key-length',type=int,choices=[13,14],default=13)
parser.add_argument('--include-refinement',action='store_true')
args=parser.parse_args()

ROOT=Path(__file__).resolve().parents[2]
OPT=ROOT/"optimization"
while not (OPT/"results/inner_full_summary.json").exists():
    time.sleep(2)
runs=json.loads((OPT/"results/inner_full_summary.json").read_text())
if args.include_refinement:
    runs+=json.loads((OPT/'results/refined_summary.json').read_text())
prefix=('refined.' if args.include_refinement else '')+f'header{args.key_length}.'
votes={};supports={}
for index,run in enumerate(runs):
    with open(run['output']) as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            body=row['payload'];bit=1<<(index*2+(row['direction']=='R'))
            votes[body]=votes.get(body,0)|bit
            supports[body]=max(supports.get(body,0),int(row['support']))
    print('loaded predictions',run['tag'],len(votes),flush=True)

# Candidate selection below uses predictions, vote counts, and seed support only.
# No reference file has been opened yet.
groups=defaultdict(list)
for body,mask in votes.items():
    groups[body[:args.key_length]].append((mask.bit_count(),body,mask))
choices=[]
for candidates in groups.values():
    candidates.sort(reverse=True)
    count,body,mask=candidates[0]
    if len(candidates)>1 and candidates[1][0]==count:
        continue
    second=candidates[1][0] if len(candidates)>1 else 0
    all_sources=0
    for _,_,source in candidates:all_sources|=source
    choices.append((body,count,all_sources.bit_count(),count-second,supports[body]))

# References now score the already-selected bodies, without repairing them.
truth=set()
with gzip.open(ROOT/'data/id20.refs.txt.gz','rt') as handle:
    for line in handle:truth.add(line.strip()[3:-3][20:130])
assert len(truth)==607150
counts=np.array([row[1] for row in choices]); totals=np.array([row[2] for row in choices])
margins=np.array([row[3] for row in choices]); support=np.array([row[4] for row in choices])
correct=np.array([row[0] in truth for row in choices])
grid=[]
for minimum_votes in range(1,int(counts.max())+1):
    for agreement in [.5,.67,.8,1.0]:
        for margin in [1,2,3]:
            gate=(counts>=minimum_votes)&(counts>=agreement*totals)&(margins>=margin)
            hist=np.bincount(support[gate],minlength=129)
            exact_hist=np.bincount(support[gate&correct],minlength=len(hist))
            cumulative=np.cumsum(hist[::-1])[::-1]
            exact_cumulative=np.cumsum(exact_hist[::-1])[::-1]
            for min_support in list(range(1,33))+[48,64,128]:
                returned=int(cumulative[min_support]);exact=int(exact_cumulative[min_support])
                grid.append({'minimum_votes':minimum_votes,'minimum_agreement':agreement,'minimum_margin':margin,
                    'minimum_seed_support':min_support,'returned':returned,'exact':exact,'wrong':returned-exact,
                    'exact_rate':exact/607150,'precision':exact/returned if returned else 0})
best=max(grid,key=lambda row:(row['exact'],-row['wrong']))
zero=[row for row in grid if row['wrong']==0 and row['returned']>0]
best_zero=max(zero,key=lambda row:row['exact']) if zero else None
precise=[row for row in grid if row['precision']>=.9999 and row['returned']>0]
balanced=max(precise,key=lambda row:(row['exact'],-row['wrong'])) if precise else None
selected={'maximum_exact':best,'precision_99_99':balanced,'zero_observed_errors':best_zero}
for name,rule in selected.items():
    if rule is None:continue
    output=OPT/'results'/f'{prefix}{name}.fa.gz'
    with gzip.open(output,'wt',compresslevel=1) as handle:
        n=0
        for body,count,total,margin,seed_support in choices:
            if count<rule['minimum_votes'] or count<rule['minimum_agreement']*total or margin<rule['minimum_margin'] or seed_support<rule['minimum_seed_support']:continue
            # IDs and sequences are predictions; no truth labels or repairs.
            handle.write(f'>pred_{body[:args.key_length]} votes={count} support={seed_support}\n{body}\n');n+=1
        assert n==rule['returned']
    rule['fasta']=str(output)
summary={'source_runs':[row['tag'] for row in runs], 'independent_sources_note':'Each parameter/direction supplies at most one vote per body. Sources are correlated, not independent evidence probabilities.',
    'candidate_selection':f'unique maximum votes per {args.key_length} nt prefix, then globally fixed vote/agreement/margin/seed-support thresholds; no per-target truth choice',
    'truth_usage':'global rule evaluation/selection only; predictions never repaired or filtered by matching a reference',
    'rules_evaluated':len(grid),'selected':selected,'rs_validation':'not_performed'}
(OPT/'results'/f'{prefix}ensemble_summary.json').write_text(json.dumps(summary,indent=2))
(OPT/'results'/f'{prefix}ensemble_rules.json').write_text(json.dumps(grid,indent=2))
print(json.dumps(summary,indent=2),flush=True)
