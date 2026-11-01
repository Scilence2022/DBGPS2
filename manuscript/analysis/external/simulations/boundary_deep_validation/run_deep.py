"""Predeclared independent validation; counts each input at both k values."""
from pathlib import Path
import csv
import hashlib
import itertools
import json
import platform
import subprocess
import time
from deep_simulate import generate

R=Path(__file__).resolve().parent
CS=[1,3,5]; RS=[1.25,1.5,2,2.55,3,4,8,0]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def jobs():
    out=[]
    def add(cohort,seed,n,geometry,L,D,q,e,channel='substitution',sampling='budget',sigma=0):
        out.append(dict(cohort=cohort,seed=seed,n=n,geometry=geometry,length=L,depth=D,
                        recovery=q,error=e,channel=channel,sampling=sampling,abundance_sigma=sigma))
    for seed in [7701,8802,9903]:
        for geo,L,D,q,e in [('regular',200,30,1,.1),('regular',45,100,.05,.01),
                            ('regular',80,100,.05,.03),('bernoulli',25,30,1,.01),
                            ('bernoulli',45,30,.05,.01),('regular',200,30,1,.03)]:
            add('holdout',seed,1024,geo,L,D,q,e)
        for D,channel in itertools.product([30,100],['substitution','burst','persistent']):
            add('error_structure',seed,1024,'regular',80,D,1,.03,channel)
        for e in [.03,.1]:
            add('mixed_holdout',seed,1024,'regular',80,100,.05,e,'mixed')
        for L in [45,80]:
            add('abundance_holdout',seed,1024,'bernoulli',L,100,.05,.01,sampling='poisson',sigma=1.5)
    for seed,n in itertools.product([17701,18802,19903],[1024,8192,65536]):
        add('scale',seed,n,'regular',45,100,.05,.01)
    for seed,geo,L,q,D in itertools.product([11011,12012,13013],['regular','bernoulli'],[25,45,80],[1,.05],[10,30,100]):
        add('poisson_theory',seed,1024,geo,L,D,q,.01,sampling='poisson')
    assert len(out)==165
    return out


def make_lock():
    files=['simulate.py','deep_simulate.py','blind_policy.py','run_deep.py','smoke_deep.py','bin/blind_profile']
    files += [str(p.relative_to(R)) for p in sorted((R/'src').iterdir()) if p.is_file()]
    return dict(design='Independent holdout, structured errors, nested scale, and Poissonized coverage theory',
                n_graphs=len(jobs())*2,k=[21,31],cutoffs=CS,ratios=RS,
                policies=['default','global_crc_r','global_crc_cr','union_consensus','two_support_consensus'],
                truth_use='Only outcome scoring and oracle diagnostics. Locked policies receive native candidates and truth-independent statuses only.',
                tie_break='Maximum accepted count; then larger cutoff for c-r policy; then smallest finite r. Off sorts last.',
                jobs=jobs(),measurement_files={p:digest(R/p) for p in files},platform=platform.platform())


def main():
    lock=make_lock(); lp=R/'design_and_source_lock.json'
    if lp.exists():assert json.loads(lp.read_text())==lock
    else:lp.write_text(json.dumps(lock,indent=2))
    assert json.loads((R/'smoke_validation.json').read_text())['passed']
    for d in ['runs','scratch','references','clean_counts']:(R/d).mkdir(exist_ok=True)
    start=time.monotonic();count=0
    for job in jobs():
        key=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:16]
        meta_path=R/'runs'/f'{key}.input.json'
        todo=[k for k in [21,31] if not (R/'runs'/f'{key}_k{k}.json').exists()]
        if todo:
            readpath=R/'scratch/reads.fa';refpath=R/'references'/f'seed{job["seed"]}_n{job["n"]}.fa'
            meta=generate(job,readpath,refpath,R/'clean_counts'/key)
            meta.update(settings=job,input_sha256=digest(readpath),reference_sha256=digest(refpath))
            if meta_path.exists():assert json.loads(meta_path.read_text())==meta
            else:meta_path.write_text(json.dumps(meta,indent=2))
        for k in [21,31]:
            prefix=R/'runs'/f'{key}_k{k}';record=Path(str(prefix)+'.json')
            if record.exists():count+=1;continue
            cmd=[str(R/'bin/blind_profile'),'-x',str(refpath),'-a','1000','-b',str(999+job['n']),
                 '-k',str(k),'-t','3','-N','-v','257' if job['n']<=1024 else '4093',
                 '-c',','.join(map(str,CS)),'-r',','.join(map(str,RS)),'-o',str(prefix),str(readpath)]
            with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as err:
                p=subprocess.run(cmd,stdout=out,stderr=err)
            assert p.returncode==0,(job,k,p.returncode)
            grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'))
            assert len(grid)==24
            for row in grid:
                assert 0<=int(row['correct'])<=int(row['ceil_adm'])<=int(row['ceil_cov'])<=job['n']
                assert int(row['search_gap'])==int(row['overflow_adm'])+int(row['ambiguous_adm'])
            record.write_text(json.dumps(dict(settings=job,k=k,command=cmd,grid=grid),indent=2))
            count+=1
            print('COMPLETE',count,'/330',job['cohort'],job['seed'],job['n'],job['geometry'],job['length'],job['depth'],job['channel'],k,'elapsed',round(time.monotonic()-start,1),flush=True)
    assert count==330
    print('DEEP_SYNTHETIC_ALL_COMPLETE',flush=True)


if __name__=='__main__':main()
