from pathlib import Path
import subprocess,os,json,sys
R=Path(__file__).resolve().parent;S=R/'smoke';S.mkdir(exist_ok=True)
sys.path.insert(0,str(R));from run_benchmark import canonicalize
paths={'nbt':('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa',20,130,29,0,1.9),'song':('/Users/song/dbgps2_ht/ref.fa',18,182,31,1,2.55),'fountain':('/Users/song/Github-Repos/DBGPS2-dnafountain/analysis/runs/fullint/dbgps_do.strands',0,152,21,0,6.)}
records=[]
for mode,(file,lo,hi,k,c,r) in paths.items():
 lines=Path(file).read_text().splitlines();seqs=[s.split('\t')[-1] for s in lines if s and not s.startswith('>')][:100];want={s[lo:hi] for s in seqs};reads=S/f'{mode}.fa';reads.write_text(''.join(f'>{j}_{i}\n{s}\n' for j in range(5) for i,s in enumerate(seqs)))
 greads=reads
 if mode=='fountain':greads=S/'fountain.graph.fa';greads.write_text(''.join(f'>{j}_{i}\nCCTGCAGAGT{s}T\n' for j in range(5) for i,s in enumerate(seqs)))
 out=S/f'{mode}.graph.raw.txt';cmd=[str(R/'bin/format_graph'),mode,str(greads),str(R/'src/legacy-2base.map'),str(k),str(c),str(r),str(out),'0'];env=os.environ.copy();env['DBGPS_DETECT_ONLY']='1';subprocess.run(cmd,check=True,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 final=S/f'{mode}.graph.txt';canonicalize(out,final,mode);got=set(final.read_text().splitlines());assert got==want,(mode,'graph',len(got),len(want),list(want-got)[:1]);print(mode,'graph 100/100',flush=True)
 prefix=S/f'{mode}.msa';subprocess.run([sys.executable,str(R/'msa_pipeline.py'),mode,str(reads),str(prefix),'--workers','2'],check=True,stdout=subprocess.DEVNULL)
 final=S/f'{mode}.msa.txt';canonicalize(Path(str(prefix)+'.valid.txt'),final,mode);got=set(final.read_text().splitlines());assert got==want,(mode,'msa',len(got),len(want));print(mode,'MSA 100/100',flush=True);records.append(dict(format=mode,designs=100,graph_exact=len(want),msa_exact=len(want)))
(R/'smoke_validation.json').write_text(json.dumps(records,indent=2))
