"""Reference-free exact-prefix clustering, bounded read sampling, MUSCLE, vote.
One MSA per eligible cluster; no reference DNA or reference-derived addresses.
"""
import heapq
import argparse,collections,concurrent.futures,json,random,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parent
TRANS=bytes.maketrans(b'ACGTN',b'TGCAN')
LEGACY={line.split()[0].encode():int(line.split()[1],2) for line in (R/'src/legacy-2base.map').read_text().splitlines() if len(line.split())==2 and line.split()[0] not in ('bases','bits') and not line.startswith('#')}
FORMATS={'nbt':(110,14,b'GCAACAAGTCAATCCGT',b'AATTGAATGCTTGCTTGCCG'), 'song':(164,16,b'GCAGAGTAGCATGTC',b'CTGACACTGATGCATCCG'),'fountain':(152,16,b'',b'')}
def bodies(s,mode):
 n,keylen,anchor,end=FORMATS[mode]
 if mode=='fountain':
  if len(s)==n and set(s)<=set(b'ACGT'):yield s,False
  return
 if mode=='nbt':
  # Normalize each molecule once, avoiding duplicate F/R MSA clusters.
  for oriented in (s,s.translate(TRANS)[::-1]):
   pos=oriented.find(anchor)
   if pos>=0:
    b=oriented[pos+len(anchor):];ep=b.find(end);b=b[:ep] if ep>=0 else b[:n]
    if keylen<=len(b)<=n+15 and set(b)<=set(b'ACGT'):yield b,False
    return
  # Rescue fully observed bodies from the reverse primer if the forward
  # anchor is damaged. The canonical prefix must itself be observed.
  for oriented in (s,s.translate(TRANS)[::-1]):
   pos=oriented.find(b'CAAGCAAGCATTCAATT')
   if pos>=0:
    b=oriented[pos+17:];ep=b.find(b'ACGGATTGACTTGTTGCACT');b=b[:ep] if ep>=0 else b[:n]
    if len(b)==n and set(b)<=set(b'ACGT'):yield b.translate(TRANS)[::-1],False
    return
  return
 anchors=[(anchor,end,False)]
 for oriented in (s,s.translate(TRANS)[::-1]):
  for a,e,reverse in anchors:
   pos=oriented.find(a)
   if pos<0:continue
   b=oriented[pos+len(a):];endpos=b.find(e)
   if endpos>=0:b=b[:endpos]
   else:
    adapter=b.find(b'AGATCGGAAGAGC')
    if adapter>=0:b=b[:adapter]
    b=b[:n]
   if keylen<=len(b)<=n+15 and set(b)<=set(b'ACGT'):
    idx=0
    for j in range(0,16,2):idx=(idx<<4)|LEGACY[b[j:j+2]]
    if 101010102<=idx<=101295684:yield b,reverse

def reconstruct(job):
 key,reads,mode=job
 f=''.join(f'>s{i}\n{s.decode()}\n' for i,s in enumerate(reads))
 t=time.perf_counter();p=subprocess.run([str(R/'bin/muscle3'),'-quiet'],input=f,capture_output=True,text=True,check=True,timeout=120,close_fds=False)
 a=[]
 for line in p.stdout.splitlines():
  if line.startswith('>'):a.append('')
  else:a[-1]+=line.strip().upper()
 assert len(a)==len(reads) and len(set(map(len,a)))==1
 v=[]
 spans=[(next((i for i,b in enumerate(seq) if b!='-'),len(seq)),len(seq.rstrip('-'))) for seq in a]
 for col_index,col in enumerate(zip(*a)):
  counts=dict.fromkeys('ACGT-N',0)
  for row_index,b in enumerate(col):
   if mode=='song' and b=='-' and (col_index<spans[row_index][0] or col_index>=spans[row_index][1]):continue
   counts[b]+=1
  if not any(counts.values()):continue
  counts['-']*=.4;b=max(counts,key=counts.get)
  if b!='-':v.append(b)
 s=''.join(v)
 if key[0]==1:s=s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
 return s if len(s)==FORMATS[mode][0] and set(s)<=set('ACGT') else None,time.perf_counter()-t

def main():
 ap=argparse.ArgumentParser();ap.add_argument('mode');ap.add_argument('reads');ap.add_argument('prefix');ap.add_argument('--workers',type=int,default=8);a=ap.parse_args()
 start=time.perf_counter();groups={};counts=collections.Counter();nr=used=0;rng=random.Random(20261019)
 with open(a.reads,'rb') as f:
  while h:=f.readline():
   s=f.readline().strip();assert h.startswith(b'>');nr+=1
   for b,reverse in bodies(s,a.mode):
    used+=1;key=(int(reverse),b[:FORMATS[a.mode][1]]);counts[key]+=1;g=groups.setdefault(key,[])
    if a.mode=='song':
     # Coverage-first sampling from observed reads, using only known format
     # length. Ties are random; designs and true addresses are never loaded.
     item=(min(len(b),164),-abs(len(b)-164),rng.random(),b)
     if len(g)<15:heapq.heappush(g,item)
     elif item>g[0]:heapq.heapreplace(g,item)
    elif len(g)<15:g.append(b)
    else:
     i=rng.randrange(counts[key])
     if i<15:g[i]=b
 eligible=sorted(k for k,c in counts.items() if c>=3)
 ct=time.perf_counter();print('CLUSTERED',nr,used,len(groups),len(eligible),ct-start,flush=True)
 calls=0;worker_elapsed=0
 with open(a.prefix+'.candidates.txt','w') as out,concurrent.futures.ProcessPoolExecutor(a.workers) as pool:
  for s,t in pool.map(reconstruct,((key,[x[-1] for x in groups[key]] if a.mode=='song' else groups[key],a.mode) for key in eligible),chunksize=64):
   calls+=1;worker_elapsed+=t
   if s:out.write(s+'\n')
   if calls%25000==0:print('MSA',calls,'/',len(eligible),round(time.perf_counter()-ct,1),flush=True)
 vt=time.perf_counter()
 cmd=[str(R/'bin/format_graph'),a.mode,'validate',str(R/'src/legacy-2base.map'),a.prefix+'.candidates.txt',a.prefix+'.valid.txt']
 subprocess.run(cmd,check=True);end=time.perf_counter()
 d=dict(reads=nr,anchored_observations=used,groups=len(groups),eligible_clusters=len(eligible),alignment_calls=calls,min_cluster=3,alignment_read_cap=15,subsets=1,seed=20261019,gap_weight=.4,workers=a.workers,cluster_seconds=ct-start,msa_vote_seconds=vt-ct,validation_seconds=end-vt,pipeline_seconds=end-start,worker_elapsed_sum=worker_elapsed,reference_access=False,read_selection='longest observed coverage up to known 164-nt region; random ties' if a.mode=='song' else 'uniform reservoir',terminal_gaps='missing observations' if a.mode=='song' else 'weighted gaps')
 Path(a.prefix+'.stages.json').write_text(json.dumps(d,indent=2));print(json.dumps(d),flush=True)
if __name__=='__main__':main()
