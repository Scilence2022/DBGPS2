#!/bin/zsh
set -e
cd "$(dirname "$0")/../data/reads"
U1="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_1.fastq.gz"
U2="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_2.fastq.gz"
sz() { stat -f%z "$1" 2>/dev/null || echo 0; }
cap_dl() { local url=$1 file=$2 cap=$3
  while [ "$(sz $file)" -lt $cap ]; do
    curl -sS --retry 10 --retry-delay 5 -C - -o $file "$url" &
    local pid=$!
    while kill -0 $pid 2>/dev/null && [ "$(sz $file)" -lt $cap ]; do sleep 15; done
    kill $pid 2>/dev/null || true; wait $pid 2>/dev/null || true
  done
  echo "$file at $(sz $file)"
}
cap_dl "$U1" ERR1816980_1.fastq.gz $((130*1024*1024))
echo "R1_AT_130MB"
cap_dl "$U2" ERR1816980_2.fastq.gz $((140*1024*1024))
echo "R2_AT_140MB"
curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_1.fastq.gz "$U1" || true
curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_2.fastq.gz "$U2" || true
echo "ALL_DONE"
