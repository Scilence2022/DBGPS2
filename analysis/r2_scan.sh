#!/bin/zsh
# Scan the local R2 against the server copy in 8 MB steps (8 KB probe each);
# print offsets where bytes differ, for targeted repair.
set -e
R2=/Users/song/Github-Repos/DBGPS2-dnafountain/data/reads/ERR1816980_2.fastq.gz
U="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR181/000/ERR1816980/ERR1816980_2.fastq.gz"
SIZE=1581555238
STEP=$((8*1024*1024))
P=8192
off=0
while [ $off -lt $SIZE ]; do
  curl -s --max-time 40 -r $off-$((off+P-1)) -o /tmp/probe.bin "$U"
  L=$(dd if=$R2 bs=1 skip=$off count=$P 2>/dev/null | shasum -a 256 | cut -d' ' -f1)
  Rr=$(shasum -a 256 /tmp/probe.bin | cut -d' ' -f1)
  if [ "$L" != "$Rr" ]; then echo "DIFF $off"; fi
  off=$((off+STEP))
done
echo SCAN_DONE
