#!/bin/zsh
# FULL dataset pipeline using pre-converted FASTA files (r1g.fa, r1pg.fa, r2g.fa).
set -e
cd "$(dirname "$0")/.."
OUT=analysis/runs/full
ls -la $OUT/r1g.fa $OUT/r1pg.fa $OUT/r2g.fa

echo "== seeds =="
python3 analysis/extract_seeds.py --r1 $OUT/r1g.fa --out $OUT/seeds_all.tsv --min_count 1
awk '{print $1}' $OUT/seeds_all.tsv > $OUT/seedlist.txt

echo "== DBGPS2 full (k21, c0, r6, detect-only, -F 33) =="
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -F 33 -I $OUT/seedlist.txt -o $OUT/dbgps -t 8 \
    $OUT/r1pg.fa $OUT/r2g.fa 2>&1 | grep -E "Strands decoded|decoding time|counting time"

echo "== paper-style full =="
python3 analysis/paper_pipeline.py --r1 $OUT/r1g.fa --r2 $OUT/r2g.fa --out $OUT/paper_droplets.txt

echo "== cross validation and union =="
python3 - <<'EOF'
import sys
sys.path.insert(0,'analysis')
from fountain_common import dna_to_bytes
paper={}
for line in open('analysis/runs/full/paper_droplets.txt'):
    s=dna_to_bytes(line.strip()[:16]); paper[int.from_bytes(s,'big')]=line.strip()
d2={}
for line in open('analysis/runs/full/dbgps.strands'):
    s,dna=line.strip().split('\t'); d2[int(s)]=dna
common=set(d2)&set(paper)
match=sum(1 for s in common if d2[s]==paper[s])
print(f"DBGPS2: {len(d2)}; paper: {len(paper)}; common {len(common)}; "
      f"exact {match} ({100*match/max(1,len(common)):.2f}%); "
      f"DBGPS2-only {len(set(d2)-set(paper))}; paper-only {len(set(paper)-set(d2))}")
seen=set()
with open('analysis/runs/full/union_droplets.txt','w') as out:
    for s,dna in d2.items():
        out.write(dna+"\n"); seen.add(s)
    extra=0
    for s,dna in paper.items():
        if s not in seen:
            out.write(dna+"\n"); extra+=1
print("union:", len(seen)+extra)
EOF

echo "== decode DBGPS2-only =="
cut -f2 $OUT/dbgps.strands > $OUT/droplets_dbgps.txt
python3 analysis/decode_fountain.py $OUT/droplets_dbgps.txt --out $OUT/decoded_dbgps.bin 2>&1 | head -2
echo "== decode union =="
python3 analysis/decode_fountain.py $OUT/union_droplets.txt --out $OUT/decoded_union.bin 2>&1 | head -2
md5 -q $OUT/decoded_union.bin $OUT/decoded_dbgps.bin
echo "expected: 8651e90d3a013178b816b63fdbb94b9b"
echo FULL_RUN_DONE
