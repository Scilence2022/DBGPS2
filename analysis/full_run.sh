#!/bin/zsh
# FULL dataset run (all ~15.9M read pairs). R2's gzip trailer check fails
# (cosmetic damage from interrupted downloads); zcat emits all data anyway.
set -e
cd "$(dirname "$0")/.."
OUT=analysis/runs/full
mkdir -p $OUT
R1=data/reads/ERR1816980_1.fastq.gz
R2=data/reads/ERR1816980_2.fastq.gz

echo "== convert to FASTA (lockstep, ACGT-only, drop pairs with N) =="
python3 - <<'EOF'
import gzip, zlib

class RawGz:
    """gzip reader that tolerates a broken trailer."""
    def __init__(self, path):
        self.f = open(path, 'rb')
        head = self.f.read(10)
        assert head[:2] == b'\x1f\x8b'
        self.d = zlib.decompressobj(-zlib.MAX_WBITS)
        self.buf = b''
    def line(self):
        while b'\n' not in self.buf:
            chunk = self.f.read(1 << 22)
            if not chunk:
                if self.buf:
                    l, self.buf = self.buf, b''
                    return l
                raise EOFError
            self.buf += self.d.decompress(chunk)
        l, self.buf = self.buf.split(b'\n', 1)
        return l

f1 = gzip.open('data/reads/ERR1816980_1.fastq.gz', 'rb')
f2 = RawGz('data/reads/ERR1816980_2.fastq.gz')
o1 = open('analysis/runs/full/r1.fa', 'w')
o1p = open('analysis/runs/full/r1p.fa', 'w')
o2 = open('analysis/runs/full/r2.fa', 'w')
P = 'CCTGCAGAGT'
n = kept = 0
try:
    while True:
        h1 = f1.readline(); h2 = f2.line()
        s1 = f1.readline(); s2 = f2.line()
        f1.readline(); f2.line()
        f1.readline(); f2.line()
        if not h1 or not h2:
            break
        n += 1
        a = s1.strip(); b = s2.strip()
        if set(a) <= set(b'ACGT') and set(b) <= set(b'ACGT'):
            name = h1.split()[0][1:].decode()
            o1.write(f'>{name}\n{a.decode()}\n')
            o1p.write(f'>{name}\n{P}{a.decode()}\n')
            o2.write(f'>{name}\n{b.decode()}\n')
            kept += 1
except EOFError:
    pass
print(f'pairs read: {n}, ACGT-clean kept: {kept}', flush=True)
EOF
wc -l $OUT/r1.fa $OUT/r2.fa

echo "== seeds =="
python3 analysis/extract_seeds.py --r1 $OUT/r1.fa --out $OUT/seeds_all.tsv --min_count 1
awk '{print $1}' $OUT/seeds_all.tsv > $OUT/seedlist.txt

echo "== DBGPS2 full (k21, c0, r6, detect-only, singleton filter -F 33) =="
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -F 33 -I $OUT/seedlist.txt -o $OUT/dbgps -t 8 \
    $OUT/r1p.fa $OUT/r2.fa 2>&1 | grep -E "Strands decoded|decoding time|counting time"

echo "== paper-style full =="
python3 analysis/paper_pipeline.py --r1 $OUT/r1.fa --r2 $OUT/r2.fa --out $OUT/paper_droplets.txt

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
