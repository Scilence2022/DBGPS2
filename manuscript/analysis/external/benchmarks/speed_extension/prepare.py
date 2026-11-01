"""Shared, reference-free paired-read merging, excluded from reconstruction timing."""
import subprocess,time,json,hashlib,gzip
from pathlib import Path
R=Path(__file__).resolve().parent
D=R/'data';D.mkdir(exist_ok=True)
datasets={'fountain':('/Users/song/Github-Repos/DBGPS2-dnafountain/data/reads/ERR1816980_1.fastq.gz','/Users/song/Github-Repos/DBGPS2-dnafountain/data/reads/ERR1816980_2.fastq.gz'),'song_ht176':('/Users/song/Documents/123/HT176_FDDP220002706-1a_1.fq.gz','/Users/song/Documents/123/HT176_FDDP220002706-1a_2.fq.gz')}
records=[]
for label,paths in datasets.items():
 dst=D/f'{label}.fa'; merged=D/f'{label}.merged.fq'
 meta=D/f'{label}.preparation.json'
 if meta.exists(): print('EXISTS',label,flush=True);continue
 start=time.perf_counter();cmd=[str(R/'bin/NGmerge'),'-1',paths[0],'-2',paths[1],'-o',str(merged),'-n','8','-d','-y','-v']
 print('MERGE',label,flush=True)
 with (D/f'{label}.merge.stderr').open('w') as err:subprocess.run(cmd,stderr=err,check=True)
 merge_elapsed=time.perf_counter()-start
 n=retained=0
 with merged.open('rb') as f,dst.open('wb') as out:
  while h:=f.readline():
   s=f.readline().strip();p=f.readline();q=f.readline().strip();assert h.startswith(b'@') and p.startswith(b'+') and len(s)==len(q)
   n+=1
   if label=='fountain' and len(s)!=152:continue
   retained+=1;out.write(b'>'+str(retained).encode()+b'\n'+s+b'\n')
 row=dict(label=label,inputs=list(paths),command=cmd,merge_seconds=merge_elapsed,preparation_seconds=time.perf_counter()-start,merged_records=n,retained_records=retained,shared_fasta=str(dst),filter='152 nt merged molecules' if label=='fountain' else 'all successfully merged molecules; arbitrary lengths',failed_merges='excluded equally for all reconstruction methods')
 meta.write_text(json.dumps(row,indent=2));records.append(row);print(json.dumps(row),flush=True)
 merged.unlink()
