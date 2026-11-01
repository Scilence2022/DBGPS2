#!/bin/zsh
# Phase A: R2 to 72MB; Phase B: R1 to 70MB; Phase C: resume R2 full; Phase D: resume R1 full.
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
cap_dl "$U2" ERR1816980_2.fastq.gz $((72*1024*1024))
echo "PHASE_A_DONE"
cap_dl "$U1" ERR1816980_1.fastq.gz $((70*1024*1024))
echo "PHASE_B_DONE"
curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_2.fastq.gz "$U2" || true
echo "R2 full: $(sz ERR1816980_2.fastq.gz)"
curl -sS --retry 20 --retry-delay 5 -C - -o ERR1816980_1.fastq.gz "$U1" || true
echo "R1 full: $(sz ERR1816980_1.fastq.gz)"
echo "ALL_DONE"
