from pathlib import Path
import subprocess,json,time
R=Path(__file__).resolve().parent
# Validate identical codec constraints on published design regions, not used in reconstruction.
for mode,path,lo,hi in [('nbt','/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa',20,130),('song','/Users/song/dbgps2_ht/ref.fa',18,182)]:
 seqs=[s.strip()[lo:hi] for s in Path(path).read_text().splitlines() if not s.startswith('>')]
 p=R/'sources'/f'{mode}_design_validation.txt';p.write_text('\n'.join(seqs)+'\n');out=Path(str(p)+'.valid')
 cmd=[str(R/'bin/format_graph'),mode,'validate',str(R/'src/legacy-2base.map'),str(p),str(out)];subprocess.run(cmd,check=True)
 accepted=out.read_text().splitlines();assert accepted==seqs,(mode,len(accepted),len(seqs));print('VALIDATED',mode,len(seqs),flush=True)
 p.unlink();out.unlink()
# Anchor-only decoration for graph ABI; payload unchanged. Known boundary lookahead T.
p=R/'data/fountain.graph.fa'
if not p.exists():
 with (R/'data/fountain.fa').open('rb') as f,p.open('wb') as out:
  while h:=f.readline():
   s=f.readline().strip();out.write(h+b'CCTGCAGAGT'+s+b'T\n')
 print('FOUNTAIN graph decoration complete',flush=True)
