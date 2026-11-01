"""Paired synthetic fragmentation/recovery channels; no reference repair of reads."""
from pathlib import Path
import binascii,hashlib,json,random
import numpy as np
P1='CCTGCAGAGTAGCATGTC';P2='CTGACACTGATGCATCCG'
def seed_of(*items):return int.from_bytes(hashlib.sha256(json.dumps(items).encode()).digest()[:8],'little')
def encode(data):return ''.join('ACGT'[(b>>s)&3] for b in data for s in [6,4,2,0])
def designs(seed,n):
 rng=random.Random(seed);out={}
 for index in range(1000,1000+n):
  raw=index.to_bytes(4,'big')+bytes(rng.randrange(256) for _ in range(35));raw+=binascii.crc_hqx(raw,0).to_bytes(2,'big')
  out[index]=P1+encode(raw)+P2
 assert all(len(s)==200 for s in out.values());return out
def fragments(ref,seed,length,depth,recovery,abundance_sigma=0):
 rows=[];physical=[];frag_lengths=[];base_depth=[]
 arng=np.random.default_rng(seed_of(seed,'abundance'));abundance=arng.lognormal(-abundance_sigma**2/2,abundance_sigma,len(ref)) if abundance_sigma else np.ones(len(ref))
 for j,(index,seq) in enumerate(ref.items()):
  rng=random.Random(seed_of(seed,length,index,'fragments'))
  target=depth*len(seq)*abundance[j];retained=0;molecules=0;diff=np.zeros(201,dtype=np.int32)
  while retained<target:
   molecules+=1
   if length>=len(seq):cuts=[0,len(seq)]
   else:
    offset=rng.randrange(1,length+1);cuts=[0,*range(offset,len(seq),length),len(seq)]
   for a,b in zip(cuts,cuts[1:]):
    if a and b<len(seq) and rng.random()>recovery:continue
    rows.append((f'{index}_{molecules}_{a}_{b}',seq[a:b]));retained+=b-a;diff[a]+=1;diff[b]-=1;frag_lengths.append(b-a)
  physical.append(molecules);base_depth.append(np.cumsum(diff)[:-1])
 depths=np.asarray(base_depth);center=depths[:,80:120].mean();ends=np.concatenate([depths[:,18:50],depths[:,150:182]],axis=1).mean()
 return rows,dict(reads=len(rows),retained_bases=sum(frag_lengths),mean_fragment_length=float(np.mean(frag_lengths)),median_fragment_length=float(np.median(frag_lengths)),mean_physical_molecules=float(np.mean(physical)),observed_base_depth=float(depths.mean()),within_strand_base_cv=float(np.mean(depths.std(axis=1)/depths.mean(axis=1))),between_strand_base_cv=float(depths.mean(axis=1).std()/depths.mean()),center_end_ratio=float(center/ends) if ends else 0)
def corrupt(rows,seed,error,channel='substitution'):
 lengths=np.array([len(s) for _,s in rows],dtype=np.int64);offsets=np.r_[0,np.cumsum(lengths)];raw=''.join(s for _,s in rows).encode();lut=np.zeros(256,dtype=np.uint8)
 for i,b in enumerate(b'ACGT'):lut[b]=i
 x=lut[np.frombuffer(raw,dtype=np.uint8)];rng=np.random.default_rng(seed);u=rng.random(len(x));extra=rng.integers(1,4,len(x),dtype=np.uint8)
 if channel=='substitution':sub=u<error;delete=np.zeros(len(x),bool);insert=delete.copy()
 else:
  assert channel=='mixed';sub=u<error*.5;delete=(u>=error*.5)&(u<error*.75);insert=(u>=error*.75)&(u<error)
 x[sub]=(x[sub]+extra[sub])%4
 repeat=1+insert.astype(np.int32)-delete.astype(np.int32);starts=np.r_[0,np.cumsum(repeat)];y=np.repeat(x,repeat);insertion_base=rng.integers(0,4,len(x),dtype=np.uint8);y[starts[:-1][insert]]=insertion_base[insert]
 encoded=np.frombuffer(b'ACGT',dtype=np.uint8)[y].tobytes();ends=starts[offsets];out=[]
 for j,(name,_) in enumerate(rows):out.append((name,encoded[int(ends[j]):int(ends[j+1])].decode()))
 return out,dict(substitutions=int(sub.sum()),deletions=int(delete.sum()),insertions=int(insert.sum()),source_bases=len(raw),output_bases=len(encoded),empty_reads=sum(not s for _,s in out))
def write_fasta(path,rows):
 with Path(path).open('w') as f:
  for name,s in rows:
   if s:f.write(f'>{name}\n{s}\n')
