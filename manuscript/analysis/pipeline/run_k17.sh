#!/bin/bash
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/ksweep
until [ -f $O/ALL_DONE ]; do sleep 20; done
for s in HT5 HT_con HT4 HT3 HT2 HT1; do
  while [ $(pgrep -f "dbgps2_profile_v2 -k" | wc -l) -ge 2 ]; do sleep 5; done
  ./dbgps2_profile_v2 -k 17 -x ref.fa -M $W/examples/legacy-2base.map \
    -c 1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100 -r 0,2,2.55,3,4 -t 8 \
    -o $O/${s}_k17 fq/${s}_1.fq fq/${s}_2.fq > $O/${s}_k17.stdout 2> $O/${s}_k17.stderr &
  sleep 2
done
wait
echo done > $O/K17_DONE
