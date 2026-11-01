#!/usr/bin/env python3
"""Independently reconcile exported union matrices, witnesses, and totals."""
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'input/antkowiak2020/results/inner_rs/parameter_union'

def main():
    summary=json.loads((OUT/'summary.json').read_text());n=summary['total_indices']
    with (OUT/'parameters.tsv').open() as handle:parameters=list(csv.DictReader(handle,delimiter='\t'))
    assert len(parameters)==summary['completed_parameters']==179
    expected_ids=[int(p['parameter_id']) for p in parameters]
    assert expected_ids==list(range(179))
    masks=[];zero_masks=[];counts=[0]*179;zero_counts=[0]*179;hist=Counter();blocks=[0]*12
    candidate_count=zero_count=rs_only=0
    with gzip.open(OUT/'candidate_success_matrix.tsv.gz','rt') as matrix,(OUT/'per_index.tsv').open() as index_table:
        assert matrix.readline().strip().split('\t')==['decoded_index']+['p'+str(p) for p in expected_ids]
        rows=csv.DictReader(index_table,delimiter='\t')
        for i,line in enumerate(matrix):
            fields=line.rstrip('\n').split('\t');assert int(fields[0])==i and len(fields)==180
            row=next(rows);assert int(row['decoded_index'])==i
            mask=zero_mask=0
            for p,code in enumerate(fields[1:]):
                assert code in ['0','1','2']
                if code!='0':mask|=1<<p;counts[p]+=1
                if code=='1':zero_mask|=1<<p;zero_counts[p]+=1
            hits=sum(code!='0' for code in fields[1:]);hist[hits]+=1
            assert hits==int(row['successful_candidate_parameters'])
            assert bool(mask)==bool(int(row['candidate_success']))
            assert bool(zero_mask)==bool(int(row['zero_candidate_success']))
            assert bool(mask and not zero_mask)==bool(int(row['rs_required_globally']))
            if mask:
                assert fields[int(row['first_candidate_parameter'])+1]!='0'
                blocks[i//16383]+=1
            masks.append(mask);zero_masks.append(zero_mask)
            candidate_count+=bool(mask);zero_count+=bool(zero_mask);rs_only+=bool(mask and not zero_mask)
        assert next(rows,None) is None and len(masks)==n
    for p,record in enumerate(parameters):
        assert counts[p]==int(record['candidate_correct'])
        assert zero_counts[p]==int(record['zero_candidate_correct'])
    assert candidate_count==summary['candidate_union_correct']
    assert zero_count==summary['zero_candidate_union_correct']
    assert rs_only==summary['requires_rs_across_all_parameters']
    assert blocks==[b['correct'] for b in summary['per_block']]
    assert {str(k):v for k,v in hist.items()}==summary['success_parameter_count_histogram']
    reconstructed=[0]*n;pairs=0;pair_zero=[0]*n
    with gzip.open(OUT/'successful_parameter_pairs.tsv.gz','rt') as handle:
        assert handle.readline().strip()=='decoded_index\tparameter_id\tcorrect_zero_paths\tcorrect_one_paths'
        for line in handle:
            i,p,z,o=map(int,line.split());assert 0<=i<n and 0<=p<179 and z>=0 and o>=0 and z+o>0
            assert not(reconstructed[i]&(1<<p)), 'Duplicate successful index/parameter pair'
            reconstructed[i]|=1<<p
            if z:pair_zero[i]|=1<<p
            pairs+=1
    assert reconstructed==masks and pair_zero==zero_masks
    policy_counts={};combined=[0]*n
    with gzip.open(OUT/'policy_success_masks.tsv.gz','rt') as handle:
        names=handle.readline().strip().split('\t')[1:];policy_counts=dict.fromkeys(names,0)
        assert sorted(names)==sorted(summary['per_policy'])
        for i,line in enumerate(handle):
            fields=line.strip().split('\t');assert int(fields[0])==i and len(fields)==len(names)+1
            for p,field in zip(names,fields[1:]):
                mask=int(field,16);assert not(mask&~masks[i])
                policy_counts[p]+=bool(mask);combined[i]|=mask
        assert i+1==n
    for p,count in policy_counts.items():assert count==summary['per_policy'][p]['correct_union']
    assert sum(bool(x) for x in combined)==summary['all_policy_output_union_correct']
    never=[int(line) for line in (OUT/'never_recovered_indices.txt').read_text().splitlines()]
    assert never==[i for i,mask in enumerate(masks) if not mask]
    result=dict(passed=True,indices=n,parameters=179,successful_index_parameter_pairs=pairs,
                candidate_union_correct=candidate_count,zero_candidate_union_correct=zero_count,
                rs_required_globally=rs_only,validated_policy_unions=policy_counts,
                assertions='Matrix, per-index table, successful-pair path counts, parameter totals, first witnesses, policy masks, missing indices, block totals, and summary agree')
    (OUT/'independent_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
