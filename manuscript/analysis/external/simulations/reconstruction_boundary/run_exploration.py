"""Locked simulation matrix and one raw-HT176 mechanism validation; resumable."""
from pathlib import Path
import csv,hashlib,itertools,json,os,subprocess,time
from simulate import *
R=Path(__file__).resolve().parent;RS=[1.25,1.5,2,2.55,3,4,8,0];CS=[0,1,3,5]
def jobs():
 out=[]
 def add(cohort,seeds,n,ls,ds,qs,es,channel='substitution',sigma=0):
  for seed,L,D,q,e in itertools.product(seeds,ls,ds,qs,es):
   if L==200 and q!=1:continue
   out.append(dict(cohort=cohort,seed=seed,n=n,length=L,depth=D,recovery=q,error=e,channel=channel,abundance_sigma=sigma))
 add('main',[1101,2202,3303],128,[25,45,80,200],[10,30,100],[1,.25,.05],[0,.01,.03,.10])
 add('independent_confirmation',[4404,5505,6606],512,[25,45,80],[30,100],[1,.05],[.01,.03])
 add('mixed_errors',[4404,5505,6606],256,[45,80,200],[30,100],[1,.05],[.03,.10],'mixed')
 add('abundance_shift',[4404,5505,6606],256,[45,80],[30,100],[1,.05],[.01,.03],sigma=.9)
 return out
def main():
 js=jobs();lock=dict(design='Fixed-width fragments from random-phase molecular partitions, shortened terminal fragments, probabilistic internal-fragment retention; all surviving fragments sequenced. Equal target retained-base budget per design, unless abundance_shift cohort.',n_synthetic_graphs=len(js)*2,k=[21,31],cutoffs=CS,ratios=RS,rows=js,interpretation='Exploratory synthetic validation, not a fitted digital twin, file-decoding result, or speed benchmark. Three independent pools per cohort; independent_confirmation and stresses are fixed before results are examined.')
 lockpath=R/'design_lock.json'
 if lockpath.exists():assert json.loads(lockpath.read_text())==lock
 else:lockpath.write_text(json.dumps(lock,indent=2))
 allrows=[];done=0;lastbase=None;base=None;ref=None;metabase=None
 (R/'runs').mkdir(exist_ok=True);(R/'scratch').mkdir(exist_ok=True);(R/'references').mkdir(exist_ok=True)
 start=time.monotonic()
 for job in js:
  key=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:16]
  basekey=tuple(job[x] for x in ['cohort','seed','n','length','depth','recovery','abundance_sigma'])
  for k in [21,31]:
   prefix=R/'runs'/f'{job["cohort"]}_{key}_k{k}';result=Path(str(prefix)+'.json')
   if result.exists():
    allrows.extend(json.loads(result.read_text())['grid']);done+=1;continue
   if basekey!=lastbase:
    ref=designs(job['seed'],job['n']);base,metabase=fragments(ref,job['seed'],job['length'],job['depth'],job['recovery'],job['abundance_sigma']);lastbase=basekey
   refpath=R/'references'/f'seed{job["seed"]}_n{job["n"]}.fa'
   if not refpath.exists():write_fasta(refpath,ref.items())
   reads,err=corrupt(base,seed_of(*basekey,job['channel'],'errors'),job['error'],job['channel']);readpath=R/'scratch/reads.fa';write_fasta(readpath,reads)
   cmd=[str(R/'bin/boundary_profile'),'-x',str(refpath),'-a','1000','-b',str(999+job['n']),'-k',str(k),'-t','3','-N','-v','257','-c',','.join(map(str,CS)),'-r',','.join(map(str,RS)),'-o',str(prefix),str(readpath)]
   t=time.monotonic()
   with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as stderr:p=subprocess.run(cmd,stdout=out,stderr=stderr)
   assert p.returncode==0,(key,k)
   grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(grid)==len(CS)*len(RS) and all(None not in z for z in grid)
   for z in grid:
    assert 0<=int(z['correct'])<=int(z['ceil_adm'])<=int(z['ceil_cov'])<=job['n']
    assert sum(int(z[f]) for f in ['lookahead_losses','seed_losses','ratio_losses','search_gap','correct'])==int(z['ceil_cov'])
    z.update(job,k=k,case_id=key,**metabase,**err)
   record=dict(settings=job,k=k,input_sha256=hashlib.sha256(readpath.read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(refpath.read_bytes()).hexdigest(),command=cmd,diagnostic_elapsed_seconds=time.monotonic()-t,grid=grid)
   result.write_text(json.dumps(record,indent=2));allrows.extend(grid);done+=1
   print('COMPLETE',done,'/',len(js)*2,job['cohort'],job['seed'],job['length'],job['depth'],job['recovery'],job['error'],k,'elapsed',round(time.monotonic()-start,1),flush=True)
   if done%10==0:write_summary(allrows)
 write_summary(allrows);print('SYNTHETIC_MATRIX_COMPLETE',done,flush=True)
 # Raw paired HT176, with no merge/trim. This is separate from the speed matrix.
 raw=[Path('/Users/song/Documents/123/HT176_FDDP220002706-1a_1.fq.gz'),Path('/Users/song/Documents/123/HT176_FDDP220002706-1a_2.fq.gz')];reference=Path('/Users/song/dbgps2_ht/ref.fa');prefix=R/'real_HT176_raw';result=Path(str(prefix)+'.json')
 if not result.exists():
  cmd=[str(R/'bin/boundary_profile'),'-x',str(reference),'-M',str(R/'src/legacy-2base.map'),'-a','101010102','-b','101295684','-k','31','-t','3','-N','-v','1009','-c','1,5','-r','1.25,1.5,2,2.55,3,4,8,0','-o',str(prefix),*map(str,raw)]
  print('REAL_HT176_START',flush=True)
  with Path(str(prefix)+'.stdout').open('w') as out,Path(str(prefix)+'.stderr').open('w') as err:p=subprocess.run(cmd,stdout=out,stderr=err)
  assert p.returncode==0
  grid=list(csv.DictReader(Path(str(prefix)+'.grid.tsv').open(),delimiter='\t'));assert len(grid)==16
  old=list(csv.DictReader(Path('/Users/song/dbgps2_ht/results/full/HT176.grid.tsv').open(),delimiter='\t'));checks=[]
  for z in grid:
   prior=next((x for x in old if float(x['r'])==float(z['r']) and x['c']==z['c']),None)
   if prior:
    for field in ['correct','wrong','no_path','overflow','no_valid','ambiguous']:assert z[field]==prior[field],(field,z,prior)
    checks.append(dict(c=z['c'],r=z['r'],matches_previous_full_raw_profile=True))
  result.write_text(json.dumps(dict(command=cmd,grid=grid,independent_existing_profile_checks=checks),indent=2))
 print('EXPLORATION_ALL_COMPLETE',flush=True)
def write_summary(rows):
 path=R/'simulation_results.csv';temp=path.with_suffix('.csv.tmp')
 with temp.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 temp.replace(path)
if __name__=='__main__':main()
