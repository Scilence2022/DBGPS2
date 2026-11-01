from pathlib import Path
import shutil,json,hashlib,subprocess
base=Path('/Users/song/.codex/worktrees/cab4/DBGPS2/output/benchmarks/2026-10-19/speed_extension');root=base.parent/'resource_control';root.mkdir(exist_ok=True);(root/'bin').mkdir(exist_ok=True)
shutil.copytree(base/'src',root/'src',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__'))
core=root/'src/fountain_DBGPS2.c';s=core.read_text();needle='typedef struct { // global data structure for kt_pipeline()';helper='''/* Benchmark-only scheduling control; reconstruction/counting rules unchanged. */
static int benchmark_pipeline_threads(void) {
    const char *s = getenv("BENCHMARK_PIPELINE_THREADS");
    int n = s ? atoi(s) : 3;
    if (n != 1 && n != 3) { fprintf(stderr, "Invalid benchmark pipeline size\\n"); exit(2); }
    return n;
}

''';assert s.count(needle)==1 and s.count('kt_pipeline(3, worker_pipeline, &pl, 3);')==2;s=s.replace(needle,helper+needle).replace('kt_pipeline(3, worker_pipeline, &pl, 3);','kt_pipeline(benchmark_pipeline_threads(), worker_pipeline, &pl, 3);');core.write_text(s)
for name in ['format_graph','format_graph_bloom']:
 p=root/'src'/f'{name}.c';s=p.read_text();needle='int main(int argc,char**argv){';assert s.count(needle)==1;s=s.replace(needle,needle+'\n const char* wenv=getenv("BENCHMARK_COUNT_WORKERS"); int benchmark_workers=wenv?atoi(wenv):8; if(benchmark_workers!=1&&benchmark_workers!=3&&benchmark_workers!=8)return 2;\n fprintf(stderr,"RESOURCE count_workers=%d pipeline_threads=%d search_workers=1\\n",benchmark_workers,benchmark_pipeline_threads());')
 s=s.replace('10000000,8,200','10000000,benchmark_workers,200');p.write_text(s)
 cmd=['cc','-O3','-o',str(root/'bin'/name),str(p),*[str(root/'src'/f) for f in ['dbgps_codec.c','dbgps_redundancy.c','dbgps_strand.c','kthread.c']],'-lz','-lpthread'];subprocess.run(cmd,check=True)
shutil.copy2(base/'bin/muscle3',root/'bin/muscle3')
p=root/'msa_pipeline.py';s=(base/'msa_pipeline.py').read_text().replace('import heapq','import contextlib,resource\nimport heapq');s=s.replace(" with open(a.prefix+'.candidates.txt','w') as out,concurrent.futures.ProcessPoolExecutor(a.workers) as pool:\n  for s,t in pool.map(reconstruct,((key,[x[-1] for x in groups[key]] if a.mode=='song' else groups[key],a.mode) for key in eligible),chunksize=64):", " with open(a.prefix+'.candidates.txt','w') as out,contextlib.ExitStack() as stack:\n  jobs=((key,[x[-1] for x in groups[key]] if a.mode=='song' else groups[key],a.mode) for key in eligible)\n  if a.workers==1: results=map(reconstruct,jobs)\n  else:\n   pool=stack.enter_context(concurrent.futures.ProcessPoolExecutor(a.workers));results=pool.map(reconstruct,jobs,chunksize=64)\n  for s,t in results:")
assert 'if a.workers==1: results=map(reconstruct,jobs)' in s;p.write_text(s)
manifest={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for folder in ['src','bin'] for p in (root/folder).rglob('*') if p.is_file()}
manifest['msa_pipeline.py']=dict(bytes=p.stat().st_size,sha256=hashlib.sha256((root/'msa_pipeline.py').read_bytes()).hexdigest())
(root/'source_manifest.json').write_text(json.dumps(manifest,indent=2));
source_script=Path(__file__).resolve()
if source_script != root/'setup_resource_control.py': shutil.copy2(source_script,root/'setup_resource_control.py')
print(root)
