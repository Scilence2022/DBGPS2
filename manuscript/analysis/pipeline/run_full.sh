#!/bin/bash
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
D=/Users/song/Documents/123
for s in HT_con HT1 HT2 HT3 HT4 HT5; do
  f1=fq/${s}_1.fq; f2=fq/${s}_2.fq
  /usr/bin/time -l ./dbgps2_profile_v2 -x ref.fa -M $W/examples/legacy-2base.map \
    -c 1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100 -r 0,1.5,2,2.55,3,4,6 -t 16 -T trace_ids.txt \
    -o results/full/$s $f1 $f2 > results/full/$s.stdout 2> results/full/$s.stderr
done
echo done > results/full/ALL_DONE
