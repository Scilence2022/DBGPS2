#!/bin/zsh
# One-tile experiment: first 375,000 read pairs (750,000 reads), matching the
# paper's headline downsampling condition. Runs DBGPS2 with several parameter
# sets, verifies, decodes, and reports.
set -e
cd "$(dirname "$0")/.."
OUT=analysis/runs/onetile
mkdir -p $OUT

echo "== convert prefixes to FASTA (375k pairs, unfiltered) =="
python3 - <<'EOF'
import gzip
def conv(path, out_path, maxn):
    n=0
    with open(out_path,'w') as out:
        try:
            with gzip.open(path,'rt') as f:
                while n<maxn:
                    h=f.readline()
                    if not h: break
                    s=f.readline().strip(); f.readline(); f.readline()
                    out.write(f">{h.split()[0][1:]}\n{s}\n"); n+=1
        except (EOFError,gzip.BadGzipFile,OSError) as e:
            print(f"{path}: truncated at {n} reads")
    return n
n1 = conv('data/reads/ERR1816980_1.fastq.gz','analysis/runs/onetile/r1.fa',375000)
n2 = conv('data/reads/ERR1816980_2.fastq.gz','analysis/runs/onetile/r2.fa',375000)
print("R1:",n1,"R2:",n2)
EOF

echo "== seeds =="
python3 analysis/extract_seeds.py --r1 $OUT/r1.fa --out $OUT/seeds_all.tsv --min_count 1
awk 'BEGIN{n=0} $2>=2{n++; print $1}' $OUT/seeds_all.tsv > $OUT/seedlist.txt
wc -l $OUT/seedlist.txt

echo "== DBGPS2 parameter sweep =="
for cfg in "A:-c 0 -r 2.55" "B:-c 0 -r 8.0 -s 2" "C:-c 1 -r 8.0 -s 2"; do
  name=${cfg%%:*}; opts=${cfg#*:}
  echo "-- config $name: $opts"
  ./DBGPS2 -k 16 -i 16 -l 128 -E rs:2B $opts -I $OUT/seedlist.txt \
      -o $OUT/dbgps_$name -t 8 $OUT/r1.fa $OUT/r2.fa 2>&1 | grep -E "Strands decoded|decoding time"
done

echo "== consensus + verification =="
python3 analysis/build_consensus.py --r1 $OUT/r1.fa --out $OUT/consensus.tsv
python3 analysis/verify_strands.py $OUT/dbgps_A.strands --consensus $OUT/consensus.tsv --report_corrected
python3 analysis/verify_strands.py $OUT/dbgps_B.strands --consensus $OUT/consensus.tsv --report_corrected
python3 analysis/verify_strands.py $OUT/dbgps_C.strands --consensus $OUT/consensus.tsv --report_corrected

echo "== fountain decode (config A) =="
cut -f2 $OUT/dbgps_A.strands > $OUT/droplets_A.txt
python3 analysis/decode_fountain.py $OUT/droplets_A.txt --out $OUT/decoded_A.bin || true
md5 $OUT/decoded_A.bin || true
echo ONETILE_DONE
