#!/usr/bin/env python3
"""Validate the globally calibrated inner-code predicate on every design."""
import gzip
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
with gzip.open(ROOT/'data/id20.refs.txt.gz','rt') as handle:
    bodies=[line.strip()[3:-3][20:130] for line in handle]
array=np.frombuffer(''.join(bodies).encode(),dtype=np.uint8).reshape(-1,110)
lookup=np.zeros(256,dtype=np.uint8)
for i,c in enumerate(b'AGCT'):lookup[c]=i
encoded=lookup[array]
digits=(encoded[:,1:]-encoded[:,:-1]-1)%4
blocks=np.sum(digits[:,:108].reshape(-1,27,4)*np.array([27,9,3,1],dtype=np.uint16),axis=2)
valid=(digits!=3).all(axis=1)&(blocks<64).all(axis=1)&(digits[:,-1]==0)
assert valid.all() and len(bodies)==607150
assert len({body[:13] for body in bodies})==607150
result={'references_checked':607150,'valid':int(valid.sum()),'invalid':int((~valid).sum()),
    'rotation_cycle':'AGCT','six_bit_blocks':27,'terminal_trit':0,'unique_13_nt_headers':607150,
    'calibration':'rotation/phase fitted globally on 512 sampled designs; no per-target repair or lookup in validator'}
(ROOT/'optimization/results/inner_validation_recheck.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
