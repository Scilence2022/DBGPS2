"""Post hoc donor attribution of false native outputs; never used for selection."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from analyze_deep import packed_reference

R=Path(__file__).resolve().parent


def unpack(word):
    return ''.join('ACGT'[(b>>s)&3] for b in word for s in [6,4,2,0])


def read_fasta(path):
    seq=[];ids=[];current=None;text=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if current is not None:ids.append(current);seq.append(''.join(text))
            current=int(line[1:].split()[0]);text=[]
        else:text.append(line)
    if current is not None:ids.append(current);seq.append(''.join(text))
    assert all(len(s)==200 for s in seq)
    return ids,seq


def main():
    f=pd.read_csv(R/'false_accept_forensics.csv')
    f=f[(f.cohort!='real_HT176_raw')&f.reference_present]
    f=f.drop_duplicates(['case_id','index','candidate_hex'])
    out=[]
    for case,part in f.groupby('case_id'):
        rec=json.loads((R/'runs'/f'{case}.json').read_text())
        ref=rec['command'][rec['command'].index('-x')+1]
        if (R/'references'/Path(ref).name).exists():ref=R/'references'/Path(ref).name
        ids,seq=read_fasta(ref);joined='\n'.join(seq);k=rec['k']
        for row in part.itertuples():
            word=unpack(bytes.fromhex(row.candidate_hex));matches={}
            anchors=list(range(0,len(word)-k+1,7))
            if anchors[-1]!=len(word)-k:anchors.append(len(word)-k)
            for q in anchors:
                token=word[q:q+k];offset=0;hits=0
                while True:
                    p=joined.find(token,offset)
                    if p<0:break
                    donor,loc=divmod(p,201);offset=p+1;hits+=1
                    if hits>5000:break
                    if loc+k>200:continue
                    left=0;right=k
                    while q+left>0 and loc+left>0 and word[q+left-1]==seq[donor][loc+left-1]:left-=1
                    while q+right<len(word) and loc+right<200 and word[q+right]==seq[donor][loc+right]:right+=1
                    key=(ids[donor],q+left,q+right,loc+left,loc+right)
                    matches[key]=right-left
            others=[(length,key) for key,length in matches.items() if key[0]!=row.index]
            best=max(others,default=(0,(None,None,None,None,None)))
            own=max([length for key,length in matches.items() if key[0]==row.index],default=0)
            out.append(dict(case_id=case,index=row.index,k=k,n=row.n,seed=row.seed,
                            candidate_hex=row.candidate_hex,nucleotide_hamming=row.nucleotide_hamming,
                            accepting_policies=row.policies_accepting_this_candidate,
                            longest_own_match=own,longest_other_match=best[0],donor_index=best[1][0],
                            candidate_match_start=best[1][1],candidate_match_end=best[1][2],
                            donor_match_start=best[1][3],donor_match_end=best[1][4]))
    pd.DataFrame(out).to_csv(R/'false_candidate_donor_traces.csv',index=False)
    print('DONOR_TRACES_COMPLETE',len(out),flush=True)


if __name__=='__main__':main()
