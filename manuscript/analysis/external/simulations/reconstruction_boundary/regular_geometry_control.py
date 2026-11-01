"""Matched regular-partition control for Bernoulli fragmentation."""
from pathlib import Path
import csv,hashlib,itertools,json,random,subprocess,time
from simulate import *
from run_exploration import RS,CS
R=Path(__file__).resolve().parent
def geometric(ref,seed,L,D,q):
 rows=[];lengths=[];physical=[];depths=[]
 for index,s in ref.items():
  rng=random.Random(seed_of(seed,L,index,'geometric'));retained=0;count=0;diff=np.zeros(201,dtype=np.int32)
  while retained<D*200:
   count+=1;cuts=[0,*[i for i in range(1,200) if rng.random()<1/L],200]
   for a,b in zip(cuts,cuts[1:]):
    if a and b<200 and rng.random()>q:continue
    rows.append((f'{index}_{count}_{a}_{b}',s[a:b]));lengths.append(b-a);retained+=b-a;diff[a]+=1;diff[b]-=1
  physical.append(count);depths.append(np.cumsum(diff)[:-1])
 depths=np.asarray(depths);ends=np.concatenate([depths[:,18:50],depths[:,150:182]],axis=1).mean()
 return rows,dict(reads=len(rows),retained_bases=sum(lengths),mean_fragment_length=float(np.mean(lengths)),median_fragment_length=float(np.median(lengths)),mean_physical_molecules=float(np.mean(physical)),observed_base_depth=float(depths.mean()),within_strand_base_cv=float(np.mean(depths.std(axis=1)/depths.mean(axis=1))),between_strand_base_cv=float(depths.mean(axis=1).std()/depths.mean()),center_end_ratio=float(depths[:,80:120].mean()/ends))
def main():
 jobs=[dict(cohort='geometry_matched_regular',seed=s,n=256,length=L,depth=D,recovery=q,error=e,channel='substitution',abundance_sigma=0) for s,L,D,q,e in itertools.product([4404,5505,6606],[25,45,80],[30,100],[1,.05],[.01,.03])]
 lock=R/'regular_geometry_design_lock.json'
 if lock.exists():assert json.loads(lock.read_text())==jobs
 else:lock.write_text(json.dumps(jobs,indent=2))
 while 'GEOMETRIC_ALL_COMPLETE' not in (R/'geometric.log').read_text()[-200:]:time.sleep(10)
 allrows=[];previous=None;done=0
 for job in jobs:
  key=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:16];basekey=tuple(job[k] for k in ['seed','length','depth','recovery'])
  for k in [21,31]:
   prefix=R/'runs'/f'regular_geometry_{key}_k{k}';recordpath=Path(str(prefix)+'.json')
   if recordpath.exists():allrows.extend(json.loads(recordpath.read_text())['grid']);done+=1;continue
   if basekey!=previous:ref=designs(job['seed'],job['n']);base,meta=fragments(ref,*basekey);previous=basekey
   refpath=R/'references'/f'seed{job["seed"]}_n256.fa'
   if not refpath.exists():write_fasta(refpath,ref.items())
   reads,errors=corrupt(base,seed_of(*basekey,'geometric_errors'),job['error']);readpath=R/'scratch/geometric.fa';write_fasta(readpath,reads)
   cmd=[str(R/'bin/boundary_profile'),'-x',str(refpath),'-a','1000','-b','1255','-k',str(k),'-t','3','-N','-v','257','-c',','.join(map(str,CS)),'-r',','.join(map(str,RS)),'-o',str(prefix),str(readpath)]
   with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as err:p=subprocess.run(cmd,stdout=out,stderr=err)
   assert p.returncode==0
   grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(grid)==32
   for z in grid:
    assert sum(int(z[x]) for x in ['lookahead_losses','seed_losses','ratio_losses','search_gap','correct'])==int(z['ceil_cov'])
    z.update(job,k=k,case_id=key,**meta,**errors)
   recordpath.write_text(json.dumps(dict(settings=job,k=k,input_sha256=hashlib.sha256(readpath.read_bytes()).hexdigest(),command=cmd,grid=grid),indent=2));allrows.extend(grid);done+=1;print('REGULAR_GEOMETRY_COMPLETE',done,'/',len(jobs)*2,flush=True)
 with (R/'regular_geometry_results.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 print('REGULAR_GEOMETRY_ALL_COMPLETE',flush=True)
if __name__=='__main__':main()
