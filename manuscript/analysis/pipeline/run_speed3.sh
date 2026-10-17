#!/bin/bash
# Quiet, serial path-search timing at c=1: DBGPS (r=1e6) versus DBGPS2 (r=2.55), no Bloom filter.
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/speed3; mkdir -p $O
for s in HT5 HT3 HT_con; do for mode in dbgps2 dbgps; do
  case $mode in dbgps) opt="-r 1000000";; dbgps2) opt="-r 2.55";; esac
  /usr/bin/time -l $W/DBGPS2 -M $W/examples/legacy-2base.map -t 3 -c 1 $opt -o $O/${s}_$mode fq/${s}_1.fq fq/${s}_2.fq \
    > $O/${s}_$mode.stdout 2> $O/${s}_$mode.time
done; done
echo done > $O/ALL_DONE
