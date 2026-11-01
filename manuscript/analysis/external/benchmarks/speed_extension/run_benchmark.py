"""Completed serial runs only. Post-timing truth scoring. Resume from result JSON."""
import argparse,hashlib,json,os,platform,re,subprocess,sys,time
from pathlib import Path
R=Path(__file__).resolve().parent
DATA={'fountain':dict(mode='fountain',n=72000,k=21,c=0,r=6.,reads=R/'data/fountain.fa',graph=R/'data/fountain.graph.fa'),
'nbt_file16':dict(mode='nbt',n=607150,k=29,c=0,r=1.9,reads=Path('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/full.fa')),
'song_ht176':dict(mode='song',n=210000,k=31,c=1,r=2.55,reads=R/'data/song_ht176.fa')}
ENV=os.environ.copy();ENV.update(DBGPS_DETECT_ONLY='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
def canonicalize(source,destination,mode):
 groups={};conflicts=set();keylen=14 if mode=='nbt' else 16
 for s in Path(source).read_text().splitlines():
  if not s:continue
  key=s[:keylen]
  if key not in groups:groups[key]=s
  elif groups[key]!=s:conflicts.add(key)
 values=[groups[k] for k in sorted(groups) if k not in conflicts]
 destination.write_text('\n'.join(values)+ ('\n' if values else ''))
 return len(conflicts)
def truth(mode):
 if mode=='fountain':return None
 p=Path('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa') if mode=='nbt' else Path('/Users/song/dbgps2_ht/ref.fa')
 lo,hi=(20,130) if mode=='nbt' else (18,182)
 return {s.strip()[lo:hi] for s in p.read_text().splitlines() if s and not s.startswith('>')}
def score(path,cfg):
 values=path.read_text().splitlines();ref=truth(cfg['mode']);exact=sum(s in ref for s in values) if ref is not None else None
 return dict(emitted=len(values),exact=exact,incorrect=len(values)-exact if exact is not None else None,recovery=exact/cfg['n'] if exact is not None else None,quality_endpoint='exact complete encoded region' if ref is not None else 'RS-zero-syndrome unique seeds; exact original design list unavailable',sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def execute(label,method,repeat,cutoff,screen=0,pilot=False):
 cfg=DATA[label];prefix=R/('pilot' if pilot else 'runs')/f'{label}_{method}_c{cutoff}_rep{repeat}';prefix.parent.mkdir(exist_ok=True)
 result=Path(str(prefix)+'.result.json')
 if result.exists():return json.loads(result.read_text())
 if method=='cluster_msa_vote':
  cmd=[sys.executable,str(R/'msa_pipeline.py'),cfg['mode'],str(cfg['reads']),str(prefix),'--workers','8'];raw=Path(str(prefix)+'.valid.txt');ratio=None
 else:
  ratio=1e6 if method=='dbgps_rule' else cfg['r'];cmd=[str(R/'bin/format_graph'),cfg['mode'],str(cfg.get('graph',cfg['reads'])),str(R/'src/legacy-2base.map'),str(cfg['k']),str(cutoff),str(ratio),str(prefix)+'.raw.txt',str(screen)];raw=Path(str(prefix)+'.raw.txt')
 Path(str(prefix)+'.load_before.txt').write_text(subprocess.check_output(['ps','-Ao','pid,pcpu,comm','-r'],text=True))
 print('START',prefix.name,flush=True);start=time.perf_counter()
 with open(str(prefix)+'.stdout','w') as out,open(str(prefix)+'.stderr','w') as err:
  p=subprocess.run(['/usr/bin/time','-l']+cmd,stdout=out,stderr=err,env=ENV)
 assert p.returncode==0,(prefix.name,p.returncode)
 dst=Path(str(prefix)+'.strands.txt');conflicts=canonicalize(raw,dst,cfg['mode']);elapsed=time.perf_counter()-start
 d=dict(dataset=label,method=method,repeat=repeat,pilot=pilot,screen_modulus=128 if screen else 1,wall_seconds=elapsed,cutoff=cutoff if ratio else None,ratio=ratio,k=cfg['k'] if ratio else None,command=cmd,exit_code=p.returncode,conflicts=conflicts,**score(dst,cfg))
 if method=='cluster_msa_vote':d['stages']=json.loads(Path(str(prefix)+'.stages.json').read_text())
 else:
  log=Path(str(prefix)+'.stderr').read_text();d['stages']={x.lower()+'_seconds':float(re.search(x+r'.*?wall=([\d.]+)',log)[1]) for x in ['COUNT','SEED','SEARCH']}
 result.write_text(json.dumps(d,indent=2));print('DONE',prefix.name,round(elapsed,2),'seconds',d['emitted'],d['exact'],d['incorrect'],flush=True);return d

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--only');ap.add_argument('--repeats',type=int,default=3);a=ap.parse_args();settings=R/'selected_settings.json'
 if a.pilot:
  chosen=json.loads(settings.read_text()) if settings.exists() else {}
  chosen.pop('song_ht5',None)
  for label in DATA:
   if a.only and label!=a.only:continue
   rows=[execute(label,'dbgps_rule',0,c,screen=1,pilot=True) for c in ([1,3,5,10,20,25,40] if label=='song_ht176' else [1,3,5,10,20])]
   best=max(rows,key=lambda r:((r['exact'] if r['exact'] is not None else r['emitted']),-(r['incorrect'] or 0),-r['cutoff']))
   chosen[label]=dict(dbgps_cutoff=best['cutoff'],criterion='highest exact count on deterministic hashed 1/128 observed-seed subset' if label!='fountain' else 'highest RS-valid count on deterministic hashed 1/128 observed-seed subset',rows=rows)
  settings.write_text(json.dumps(chosen,indent=2));return
 chosen=json.loads(settings.read_text());records=[]
 metadata=dict(started=time.strftime('%Y-%m-%dT%H:%M:%S%z'),platform=platform.platform(),python=sys.version,workers=8,repeats=a.repeats,order='serial; cyclic method rotation per dataset and repeat',timing='process start through reconstruction, validation, conflict rejection and sorted output writing; excludes common preprocessing, tuning, scoring, outer decoding',graph_baseline='absolute-cutoff DBGPS rule in same adapted core; not original published binary',bloom=False,preprocessing='NGmerge v0.3 for Fountain and Song; full untrimmed FASTA conversion for nbt; all methods share retained molecules',msa='one MUSCLE 3.8.31 alignment per exact-prefix cluster; minimum 3 observations; deterministic reservoir cap 15; weighted gap vote 0.4; format inner check; not published LSH/BMA pipeline',cpu_affinity=False,dedicated_host=False)
 (R/'benchmark_metadata.json').write_text(json.dumps(metadata,indent=2))
 for rep in range(1,a.repeats+1):
  for di,label in enumerate(DATA):
   methods=['dbgps_rule','dbgps2','cluster_msa_vote'];shift=(rep-1+di)%3;methods=methods[shift:]+methods[:shift]
   for m in methods:
    c=chosen[label]['dbgps_cutoff'] if m=='dbgps_rule' else DATA[label]['c'];records.append(execute(label,m,rep,c));(R/'results.json').write_text(json.dumps(records,indent=2))
 print('ALL_DONE',flush=True)
if __name__=='__main__':main()
