#!/bin/bash
# Quiet, serial timing of the query Bloom filter: DBGPS (r=1e6) and DBGPS2 (r=2.55), with/without -B 12.
# One counting run per configuration; three consecutive cutoffs (c=10,11,12) serve as timing replicates.
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/speed2; mkdir -p $O
for s in HT5 HT3 HT_con; do for mode in dbgps2 dbgps2_B12 dbgps dbgps_B12; do
  case $mode in
    dbgps) opt="-r 1000000";; dbgps_B12) opt="-r 1000000 -B 12";;
    dbgps2) opt="-r 2.55";; dbgps2_B12) opt="-r 2.55 -B 12";;
  esac
  /usr/bin/time -l $W/DBGPS2 -M $W/examples/legacy-2base.map -t 3 -c 10 -C 12 $opt -o $O/${s}_$mode fq/${s}_1.fq fq/${s}_2.fq \
    > $O/${s}_$mode.stdout 2> $O/${s}_$mode.time
done; done
echo done > $O/ALL_DONE
