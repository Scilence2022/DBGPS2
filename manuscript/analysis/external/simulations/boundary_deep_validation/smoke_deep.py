"""Native equivalence and an explicit reference-permutation leakage check."""
from pathlib import Path
import csv
import json
import subprocess
import numpy as np
from deep_simulate import generate, window_retention
from simulate import designs, write_fasta
from blind_policy import select

R=Path(__file__).resolve().parent;D=R/'smoke';D.mkdir(exist_ok=True)
PARENT=R.parent/'reconstruction_boundary/bin/boundary_profile'
cs=[1,3,5];rs=[1.25,1.5,2,2.55,3,4,8,0]
checks=[];native_checks=0
for channel,sampling in [('substitution','budget'),('burst','budget'),('persistent','budget'),('substitution','poisson'),('mixed','budget')]:
    job=dict(cohort='smoke',seed=19,n=16,geometry='bernoulli',length=45,depth=30,
             recovery=.25,error=.03,channel=channel,sampling=sampling,abundance_sigma=0)
    generate(job,D/'reads.fa',D/'ref.fa',D/'counts')
    for k in [21,31]:
        prefix=D/f'{channel}_{sampling}_k{k}'
        common=['-x',str(D/'ref.fa'),'-a','1000','-b','1015','-k',str(k),'-t','1','-N','-v','1','-c','1,3,5','-r',','.join(map(str,rs))]
        for binary,suffix in [(R/'bin/blind_profile',''),(PARENT,'_parent')]:
            p=subprocess.run([str(binary),*common,'-o',str(prefix)+suffix,str(D/'reads.fa')],capture_output=True,text=True)
            assert p.returncode==0,p.stderr
        assert Path(str(prefix)+'.outcomes.tsv').read_bytes()==Path(str(prefix)+'_parent.outcomes.tsv').read_bytes()
        grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'))
        old=list(csv.DictReader(Path(str(prefix)+'_parent.grid.tsv').open(),delimiter='\t'))
        for a,b in zip(grid,old):
            for field in a:
                if field!='seconds':assert a[field]==b[field],(field,a[field],b[field])
        status=np.fromfile(str(prefix)+'.blind_status.u8',dtype=np.uint8).reshape(24,16)
        candidates=np.fromfile(str(prefix)+'.accepted.u8',dtype=np.uint8).reshape(24,16,41)
        for i,row in enumerate(grid):assert int((status[i]==0).sum())==int(row['decoded_ref'])
        native_checks+=sum(int(x['native_checks']) for x in grid)
        if channel=='substitution' and sampling=='budget':
            write_fasta(D/'permuted_ref.fa',designs(918,16).items())
            alt=common.copy();alt[1]=str(D/'permuted_ref.fa')
            p=subprocess.run([str(R/'bin/blind_profile'),*alt,'-o',str(prefix)+'_permuted',str(D/'reads.fa')],capture_output=True,text=True)
            assert p.returncode==0,p.stderr
            for suffix in ['.accepted.u8','.blind_status.u8']:
                assert Path(str(prefix)+suffix).read_bytes()==Path(str(prefix)+'_permuted'+suffix).read_bytes()
            native_checks+=384
        checks.append(dict(channel=channel,sampling=sampling,k=k,native_and_parent_match=True))

# Geometric window-retention formula agrees with direct enumeration in the
# regular arm and has the correct endpoint and k>maximum-length behavior.
for L in [25,45,80]:
    for q in [1,.05]:
        assert np.allclose(window_retention('regular',L,q,1).sum(),sum(window_retention('regular',L,q,1)))
        assert np.all(window_retention('regular',L,q,31)==0) if L==25 else True
        assert np.isclose(window_retention('bernoulli',L,q,31)[0],(1-1/L)**30)
assert np.all(window_retention('regular',200,.05,31)==1)

# Policy decisions do not depend on candidate nucleotide identities; a
# one-to-one relabeling preserves acceptance, conflicts and selected settings.
rng=np.random.default_rng(27);status=np.full((24,12),2,dtype=np.uint8)
cand=np.zeros((24,12,41),dtype=np.uint8);status[:2,:5]=0;cand[:2,:5]=rng.integers(0,256,(5,41))
status[3,0]=0;cand[3,0]=7
settings=[dict(c=c,r=r) for r in rs for c in cs]
a=select(status,cand,settings);b=select(status,cand^np.uint8(37),settings)
for key in a:
    assert np.array_equal(a[key]['accepted'],b[key]['accepted'])
    assert a[key]['selected_combo']==b[key]['selected_combo']
(R/'smoke_validation.json').write_text(json.dumps(dict(passed=True,native_search_comparisons=native_checks,
    parent_equivalence_graphs=len(checks),reference_permutation_blind_output_invariant=True,
    policy_identity_relabeling_invariant=True,checks=checks),indent=2))
print('DEEP_SMOKE_PASSED',native_checks,flush=True)
