import hashlib,json,sys
from pathlib import Path
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R));from run_benchmark import DATA
paths=set()
for cfg in DATA.values():paths.add(cfg['reads']);paths.add(cfg.get('graph',cfg['reads']))
for label in ['fountain','song_ht176']:
 d=json.loads((R/f'data/{label}.preparation.json').read_text());paths.update(map(Path,d['inputs']))
paths.update([Path('/Users/song/dbgps2_ht/ref.fa'),Path('/Users/song/.codex/worktrees/bd05/DBGPS2/experiments/nbt4079/data/references.fa')])
records=[]
for p in sorted(paths):
 sh=hashlib.sha256();md=hashlib.md5()
 with p.open('rb') as f:
  while b:=f.read(8*1024*1024):sh.update(b);md.update(b)
 row=dict(path=str(p),bytes=p.stat().st_size,sha256=sh.hexdigest(),md5=md.hexdigest());records.append(row);print(p.name,row['md5'],flush=True)

(R/'input_provenance.json').write_text(json.dumps(records,indent=2))
