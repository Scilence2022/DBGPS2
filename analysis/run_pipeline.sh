#!/bin/zsh
# End-to-end DBGPS2 fragment-recovery pipeline for the DNA Fountain dataset.
# Usage: run_pipeline.sh <tag> <r1_prefix.fa> <r2_prefix.fa> [min_count]
# Requires R2 reads already reverse-complemented into <r2_prefix.fa>? No:
# DBGPS2 canonicalizes k-mers, so raw R2 FASTA works directly.
set -e
cd "$(dirname "$0")/.."
TAG=$1; R1=$2; R2=$3; MINC=${4:-2}
OUT=analysis/runs/$TAG
mkdir -p $OUT

echo "== [$TAG] seed extraction =="
python3 analysis/extract_seeds.py --r1 $R1 --out $OUT/seeds.tsv --min_count $MINC
cut -f1 $OUT/seeds.tsv > $OUT/seedlist.txt

echo "== [$TAG] DBGPS2 k-mer counting + strand search (rs:2B) =="
./DBGPS2 -k 16 -i 16 -l 128 -E rs:2B -I $OUT/seedlist.txt \
    -o $OUT/dbgps -t 8 $R1 $R2 > $OUT/dbgps_stdout.log 2>&1
tail -4 $OUT/dbgps.log

echo "== [$TAG] verification =="
python3 analysis/verify_strands.py $OUT/dbgps.strands

echo "== [$TAG] consensus =="
python3 analysis/build_consensus.py --r1 $R1 --out $OUT/consensus.tsv
python3 - <<EOF
import sys
sys.path.insert(0,'analysis')
cons={}
for line in open('$OUT/consensus.tsv'):
    seed,seq,d,m=line.split('\t'); cons[int(seed)]=seq
n=exact=mm=0
for line in open('$OUT/dbgps.strands'):
    seed,dna=line.strip().split('\t')
    if int(seed) in cons:
        n+=1
        d=sum(a!=b for a,b in zip(dna, cons[int(seed)][:152]))
        exact+= (d==0); mm+= (d>0)
print(f"vs consensus: compared {n}, exact {exact} ({100*exact/max(n,1):.2f}%), mismatched {mm}")
EOF

echo "== [$TAG] fountain decoding =="
cut -f2 $OUT/dbgps.strands > $OUT/droplets.txt
python3 analysis/decode_fountain.py $OUT/droplets.txt --out $OUT/decoded.bin
md5 $OUT/decoded.bin
