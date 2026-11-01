#!/bin/bash
cd /Users/song/dbgps2_ht
O=results/dstoreqc; mkdir -p $O
for r in 0 2.55; do
  /Users/song/Github-Repos/DBGPS-analyzer/dstoreqc analyze -k 31 -t 4 -c 0 -C 100 -r $r -o $O/HT176_r$r ref_core.fa fq/HT176_1.fq fq/HT176_2.fq > $O/HT176_r$r.stdout 2> $O/HT176_r$r.stderr
done
echo done > $O/HT176_DONE
