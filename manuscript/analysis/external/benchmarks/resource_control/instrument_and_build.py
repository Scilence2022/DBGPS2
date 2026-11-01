from pathlib import Path
import subprocess,json,hashlib
R=Path(__file__).resolve().parent
for stem in ['format_graph','format_graph_bloom']:
 p=R/'src'/f'{stem}.c';s=p.read_text()
 if '#include <sys/resource.h>' not in s:
  s=s.replace('#include <sys/time.h>','#include <sys/time.h>\n#include <sys/resource.h>\nstatic double cpu_seconds(void){struct rusage r;getrusage(RUSAGE_SELF,&r);return r.ru_utime.tv_sec+r.ru_utime.tv_usec/1e6+r.ru_stime.tv_sec+r.ru_stime.tv_usec/1e6;}')
  s=s.replace('double start=wall();','double start=wall(),count_cpu_start=cpu_seconds();').replace('"COUNT wall=%.6f\\n",wall()-start','"COUNT wall=%.6f cpu=%.6f\\n",wall()-start,cpu_seconds()-count_cpu_start')
  s=s.replace('double seedstart=wall();','double seedstart=wall(),seed_cpu_start=cpu_seconds();').replace('"SEED n=%zu wall=%.6f\\n",ns,wall()-seedstart','"SEED n=%zu wall=%.6f cpu=%.6f\\n",ns,wall()-seedstart,cpu_seconds()-seed_cpu_start')
  s=s.replace('double searchstart=wall();','double searchstart=wall(),search_cpu_start=cpu_seconds();').replace('"SEARCH queries=%lu returned=%lu wall=%.6f\\n",queries,returned,wall()-searchstart','"SEARCH queries=%lu returned=%lu wall=%.6f cpu=%.6f\\n",queries,returned,wall()-searchstart,cpu_seconds()-search_cpu_start')
  if stem.endswith('bloom'):
   s=s.replace('double buildstart=wall();','double buildstart=wall(),bloom_cpu_start=cpu_seconds();').replace('filter_shift=%d wall=%.6f\\n",query_bits,h->bloom?h->bloom->n_bytes:0,filter_shift,wall()-buildstart','filter_shift=%d wall=%.6f cpu=%.6f\\n",query_bits,h->bloom?h->bloom->n_bytes:0,filter_shift,wall()-buildstart,cpu_seconds()-bloom_cpu_start')
  assert s.count('cpu=%.6f')==(4 if stem.endswith('bloom') else 3)
  p.write_text(s)
 cmd=['cc','-O3','-o',str(R/'bin'/stem),str(p),*[str(R/'src'/f) for f in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']],'-lz','-lpthread'];subprocess.run(cmd,check=True)
p=R/'msa_pipeline.py';s=p.read_text()
if 'def cpu_usage()' not in s:
 s=s.replace('def main():','def cpu_usage():\n a=resource.getrusage(resource.RUSAGE_SELF);b=resource.getrusage(resource.RUSAGE_CHILDREN)\n return a.ru_utime+a.ru_stime+b.ru_utime+b.ru_stime\n\ndef main():')
 s=s.replace('start=time.perf_counter();groups={}','start=time.perf_counter();cpu_start=cpu_usage();groups={}')
 s=s.replace('ct=time.perf_counter();print','ct=time.perf_counter();cluster_cpu_end=cpu_usage();print')
 s=s.replace('vt=time.perf_counter()','vt=time.perf_counter();msa_cpu_end=cpu_usage()')
 s=s.replace('subprocess.run(cmd,check=True);end=time.perf_counter()','subprocess.run(cmd,check=True);end=time.perf_counter();validation_cpu_end=cpu_usage()')
 s=s.replace('d=dict(reads=nr,','d=dict(cluster_cpu_seconds=cluster_cpu_end-cpu_start,msa_vote_cpu_seconds=msa_cpu_end-cluster_cpu_end,validation_cpu_seconds=validation_cpu_end-msa_cpu_end,pipeline_cpu_seconds=validation_cpu_end-cpu_start,reads=nr,')
 p.write_text(s)
manifest={str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for folder in ['src','bin'] for p in (R/folder).rglob('*') if p.is_file()}
manifest['msa_pipeline.py']=dict(bytes=(R/'msa_pipeline.py').stat().st_size,sha256=hashlib.sha256((R/'msa_pipeline.py').read_bytes()).hexdigest());(R/'source_manifest.json').write_text(json.dumps(manifest,indent=2));print('RESOURCE_TOOLS_BUILT')
