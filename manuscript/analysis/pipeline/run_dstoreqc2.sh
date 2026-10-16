#!/bin/bash
cd /Users/song/dbgps2_ht
O=results/dstoreqc; mkdir -p $O
for s in HT_con HT1 HT2 HT3 HT4 HT5; do
  f1=fq/${s}_1.fq; f2=fq/${s}_2.fq
  for r in 0 2.55; do
    /usr/bin/time -l /Users/song/Github-Repos/DBGPS-analyzer/dstoreqc analyze -k 31 -t 16 -c 0 -C 100 -r $r \
      -o $O/${s}_r$r ref_core.fa $f1 $f2 > $O/${s}_r$r.stdout 2> $O/${s}_r$r.stderr
  done
done
echo done > $O/ALL_DONE
