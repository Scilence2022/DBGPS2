#!/bin/zsh
# Sole owner of all downloading (strictly one curl at a time):
# phase 1: grow R1 prefix to R1_TARGET
# phase 2: grow R2 prefix to R2_TARGET
# phase 3: resume R1 uncapped (for later, larger experiments)
set -e
cd "$(dirname "$0")/../data/reads"
R1_TARGET=$((33*1024*1024))
R2_TARGET=$((37*1024*1024))
U1="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_1.fastq.gz"
U2="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_2.fastq.gz"

size1() { [ -f ERR1816980_1.fastq.gz ] && stat -f%z ERR1816980_1.fastq.gz || echo 0; }
size2() { [ -f ERR1816980_2.fastq.gz ] && stat -f%z ERR1816980_2.fastq.gz || echo 0; }

dl_to_cap() { # url file cap
  local url=$1 file=$2 cap=$3
  while [ "$(stat -f%z $file 2>/dev/null || echo 0)" -lt $cap ]; do
    curl -sS --retry 10 --retry-delay 5 -C - -o $file "$url" &
    local pid=$!
    while kill -0 $pid 2>/dev/null && [ "$(stat -f%z $file 2>/dev/null || echo 0)" -lt $cap ]; do
      sleep 15
    done
    kill $pid 2>/dev/null || true
    wait $pid 2>/dev/null || true
  done
  echo "$file reached $(stat -f%z $file) bytes"
}

dl_to_cap "$U1" ERR1816980_1.fastq.gz $R1_TARGET
dl_to_cap "$U2" ERR1816980_2.fastq.gz $R2_TARGET
echo "BOTH PREFIXES READY"

curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_1.fastq.gz "$U1" || true
echo "R1 full finished at $(size1) bytes"
curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_2.fastq.gz "$U2" || true
echo "R2 full finished at $(size2) bytes"
echo "ALL_DOWNLOADS_DONE"
