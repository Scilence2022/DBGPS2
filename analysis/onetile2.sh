#!/bin/zsh
# Paper one-tile condition: 750,000 read PAIRS (1.5M reads).
set -e
cd "$(dirname "$0")/.."
OUT=analysis/runs/tile750k
mkdir -p $OUT
PAIRS=750000
PREFIX10=CCTGCAGAGT

echo "== convert =="
python3 - $PAIRS <<'EOF'
import gzip, sys
maxn=int(sys.argv[1])
def conv(path, out_path, maxn, prefix=''):
    n=0
    with open(out_path,'w') as out:
        try:
            with gzip.open(path,'rt') as f:
                while n<maxn:
                    h=f.readline()
                    if not h: break
                    s=f.readline().strip(); f.readline(); f.readline()
                    out.write(f">{h.split()[0][1:]}\n{prefix}{s}\n"); n+=1
        except (EOFError,gzip.BadGzipFile,OSError) as e:
            print(f"{path}: truncated at {n}")
    return n
n1 = conv('data/reads/ERR1816980_1.fastq.gz','analysis/runs/tile750k/r1p.fa',maxn,'CCTGCAGAGT')
n2 = conv('data/reads/ERR1816980_2.fastq.gz','analysis/runs/tile750k/r2.fa',maxn)
print("R1:",n1,"R2:",n2)
# plain R1 for the paper pipeline
n=0
with gzip.open('data/reads/ERR1816980_1.fastq.gz','rt') as f, open('analysis/runs/tile750k/r1.fa','w') as out:
    try:
        while n<maxn:
            h=f.readline()
            if not h: break
            s=f.readline().strip(); f.readline(); f.readline()
            out.write(f">{h.split()[0][1:]}\n{s}\n"); n+=1
    except (EOFError,gzip.BadGzipFile,OSError): pass
print("R1 plain:",n)
EOF

echo "== seeds =="
python3 analysis/extract_seeds.py --r1 $OUT/r1.fa --out $OUT/seeds_all.tsv --min_count 1
awk '{print $1}' $OUT/seeds_all.tsv > $OUT/seedlist.txt

echo "== DBGPS2 (k21 r6 detect-only) =="
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -I $OUT/seedlist.txt -o $OUT/dbgps -t 8 $OUT/r1p.fa $OUT/r2.fa 2>&1 | grep -E "Strands decoded|decoding time"

echo "== paper-style =="
python3 analysis/paper_pipeline.py --r1 $OUT/r1.fa --r2 $OUT/r2.fa --out $OUT/paper_droplets.txt

echo "== cross validation =="
python3 - <<'EOF'
import sys
sys.path.insert(0,'analysis')
from fountain_common import dna_to_bytes
paper={}
for line in open('analysis/runs/tile750k/paper_droplets.txt'):
    s=dna_to_bytes(line.strip()[:16]); paper[int.from_bytes(s,'big')]=line.strip()
d2={}
for line in open('analysis/runs/tile750k/dbgps.strands'):
    s,dna=line.strip().split('\t'); d2[int(s)]=dna
common=set(d2)&set(paper)
match=sum(1 for s in common if d2[s]==paper[s])
print(f"DBGPS2: {len(d2)}; paper: {len(paper)}; common {len(common)}; exact {match} "
      f"({100*match/max(1,len(common)):.2f}%); DBGPS2-only {len(set(d2)-set(paper))}; "
      f"paper-only {len(set(paper)-set(d2))}")
# union for decoding, DBGPS2 preferred
seen=set()
with open('analysis/runs/tile750k/union_droplets.txt','w') as out:
    for s,dna in d2.items():
        out.write(dna+"\n"); seen.add(s)
    for s,dna in paper.items():
        if s not in seen: out.write(dna+"\n")
print("union droplets:", len(seen)+len([s for s in paper if s not in seen]))
EOF

echo "== decode: DBGPS2 only =="
cut -f2 $OUT/dbgps.strands > $OUT/droplets_dbgps.txt
python3 analysis/decode_fountain.py $OUT/droplets_dbgps.txt --out $OUT/decoded_dbgps.bin 2>&1 | head -2
echo "== decode: union =="
python3 analysis/decode_fountain.py $OUT/union_droplets.txt --out $OUT/decoded_union.bin 2>&1 | head -2
md5 -q $OUT/decoded_union.bin; echo "expected: 8651e90d3a013178b816b63fdbb94b9b"
echo TILE750K_DONE
