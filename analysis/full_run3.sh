#!/bin/zsh
# Definitive full-data run: both mates complete and md5-verified.
set -e
R=/Users/song/Github-Repos/DBGPS2-dnafountain
cd $R
OUT=analysis/runs/fullint
mkdir -p $OUT

echo "== convert (R1 plain / R1 prefixed / R2 plain / R2 'A'-prepended) =="
gzip -dc $R/data/reads/ERR1816980_1.fastq.gz 2>/dev/null | awk 'NR%4==1{print ">" NR} NR%4==2{print}' > /tmp/f_r1.fa
gzip -dc $R/data/reads/ERR1816980_1.fastq.gz 2>/dev/null | awk 'NR%4==1{print ">" NR} NR%4==2{print "CCTGCAGAGT" $0}' > /tmp/f_r1p.fa
gzip -dc $R/data/reads/ERR1816980_2.fastq.gz 2>/dev/null | awk 'NR%4==1{print ">" NR} NR%4==2{print}' > /tmp/f_r2.fa
gzip -dc $R/data/reads/ERR1816980_2.fastq.gz 2>/dev/null | awk 'NR%4==1{print ">" NR} NR%4==2{print "A" $0}' > /tmp/f_r2A.fa
awk 'END{print "R1 records:", NR/2}' /tmp/f_r1.fa
awk 'END{print "R2 records:", NR/2}' /tmp/f_r2.fa

echo "== seeds =="
python3 analysis/extract_seeds.py --r1 /tmp/f_r1.fa --out $OUT/seeds_all.tsv --min_count 1
awk '{print $1}' $OUT/seeds_all.tsv > $OUT/seedlist.txt

echo "== DBGPS2 detect-only =="
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -I $OUT/seedlist.txt -o $OUT/dbgps_do -t 8 /tmp/f_r1p.fa /tmp/f_r2A.fa 2>&1 | grep -E "Strands decoded|decoding time"

echo "== paper-style =="
python3 analysis/paper_pipeline.py --r1 /tmp/f_r1.fa --r2 /tmp/f_r2.fa --out $OUT/paper_droplets.txt

echo "== cross validation =="
python3 - <<'EOF'
import sys
sys.path.insert(0,'analysis')
from fountain_common import dna_to_bytes
paper={}
for line in open('analysis/runs/fullint/paper_droplets.txt'):
    s=dna_to_bytes(line.strip()[:16]); paper[int.from_bytes(s,'big')]=line.strip()
d2={}
for line in open('analysis/runs/fullint/dbgps_do.strands'):
    s,dna=line.strip().split('\t'); d2[int(s)]=dna
common=set(d2)&set(paper)
match=sum(1 for s in common if d2[s]==paper[s])
print(f"DBGPS2: {len(d2)}; paper: {len(paper)}; common {len(common)}; "
      f"exact {match} ({100*match/max(1,len(common)):.3f}%); "
      f"DBGPS2-only {len(set(d2)-set(paper))}; paper-only {len(set(paper)-set(d2))}")
EOF

echo "== decode DBGPS2-only (no screen) =="
cut -f2 $OUT/dbgps_do.strands > $OUT/droplets_do.txt
python3 analysis/decode_fountain.py $OUT/droplets_do.txt --no_screen --out $OUT/decoded_do.bin 2>&1 | head -2
md5 -q $OUT/decoded_do.bin
echo "expected: 8651e90d3a013178b816b63fdbb94b9b"
echo FULLINT_DONE
