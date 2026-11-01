"""New independent pools for the post hoc downcrossing-bound derivation."""
from pathlib import Path
import itertools,json,hashlib,subprocess,csv
from deep_simulate import generate
from run_deep import digest,CS,RS

R=Path(__file__).resolve().parent
jobs=[]
for seed in [21021,22022,23023]:
    for geo,L,D,q,e in [('bernoulli',80,10,.05,.01),('bernoulli',45,10,1,.03),
                        ('bernoulli',25,30,1,.01),('regular',45,30,.05,.01),
                        ('regular',45,100,.05,.01),('bernoulli',45,30,.05,.03)]:
        jobs.append(dict(cohort='crossing_confirmation',seed=seed,n=4096,geometry=geo,length=L,depth=D,
                         recovery=q,error=e,channel='substitution',sampling='poisson',abundance_sigma=0))
lock=dict(jobs=jobs,k=[21,31],cutoffs=CS,ratios=RS,
          interpretation='Prospective new-pool check of a downcrossing bound derived after inspecting the original independent-window failure.',
          files={str(p.relative_to(R)):digest(p) for p in [Path(__file__),R/'crossing_bound.py',R/'deep_simulate.py',R/'bin/blind_profile']})
lp=R/'crossing_design_and_source_lock.json'
if lp.exists():assert json.loads(lp.read_text())==lock
else:lp.write_text(json.dumps(lock,indent=2))
(R/'crossing_runs').mkdir(exist_ok=True)
for job in jobs:
    key=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:16]
    if all((R/'crossing_runs'/f'{key}_k{k}.json').exists() for k in [21,31]):continue
    ref=R/'references'/f'seed{job["seed"]}_n4096.fa';reads=R/'scratch/crossing.fa'
    meta=generate(job,reads,ref,R/'clean_counts'/key)
    meta.update(settings=job,input_sha256=digest(reads),reference_sha256=digest(ref))
    (R/'crossing_runs'/f'{key}.input.json').write_text(json.dumps(meta,indent=2))
    for k in [21,31]:
        prefix=R/'crossing_runs'/f'{key}_k{k}'
        cmd=[str(R/'bin/blind_profile'),'-x',str(ref),'-a','1000','-b','5095','-k',str(k),'-t','3','-N','-v','1021',
             '-c',','.join(map(str,CS)),'-r',','.join(map(str,RS)),'-o',str(prefix),str(reads)]
        with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as err:
            p=subprocess.run(cmd,stdout=out,stderr=err)
        assert p.returncode==0,(job,k)
        grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(grid)==24
        Path(str(prefix)+'.json').write_text(json.dumps(dict(settings=job,k=k,command=cmd,grid=grid),indent=2))
        print('CROSSING_COMPLETE',seed,job['geometry'],job['length'],job['depth'],job['error'],k,flush=True)
print('CROSSING_VALIDATION_ALL_COMPLETE',flush=True)
