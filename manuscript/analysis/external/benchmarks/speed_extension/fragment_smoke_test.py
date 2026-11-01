"""Test missing-terminal-gap handling without using references in reconstruction."""
from pathlib import Path
import subprocess,sys,json
R=Path(__file__).resolve().parent;S=R/'smoke';S.mkdir(exist_ok=True)
seqs=[s for s in Path('/Users/song/dbgps2_ht/ref.fa').read_text().splitlines() if not s.startswith('>')][:100]
p=S/'fragmented_song.fa';p.write_text(''.join(f'>{i}_{j}\n{s if j==0 else s[:50+(j%8)*10]}\n' for i,s in enumerate(seqs) for j in range(20)))
prefix=S/'fragmented_song';subprocess.run([sys.executable,str(R/'msa_pipeline.py'),'song',str(p),str(prefix),'--workers','2'],check=True,stdout=subprocess.DEVNULL)
got=set(Path(str(prefix)+'.valid.txt').read_text().splitlines());want={s[18:182] for s in seqs};assert got==want,(len(got),len(want))
(R/'fragment_smoke_validation.json').write_text(json.dumps(dict(designs=100,full_observations_per_design=1,short_prefix_fragments_per_design=19,exact=100,terminal_gaps='missing observations',true_sequences_accessed_by_reconstruction=False),indent=2));print('Fragmented MSA smoke test: 100/100 exact',flush=True)
