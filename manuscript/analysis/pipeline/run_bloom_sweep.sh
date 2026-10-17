#!/bin/bash
# Parameter sweep of DBGPS2's two Bloom modes (DBGPS has none): -B bits/key (query speed) and -F log2 bitmap bits (memory).
# DBGPS2 with r=2.55, c=1, 3 threads; one run at a time.
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/bloom_sweep; mkdir -p $O
for s in HT5 HT_con; do for mode in B8 B16 B24 F32 F36; do
  case $mode in B*) opt="-B ${mode#B}";; F*) opt="-F ${mode#F}";; esac
  /usr/bin/time -l $W/DBGPS2 -M $W/examples/legacy-2base.map -t 3 -c 1 -r 2.55 $opt -o $O/${s}_$mode fq/${s}_1.fq fq/${s}_2.fq \
    > $O/${s}_$mode.stdout 2> $O/${s}_$mode.time
done; done
echo done > $O/ALL_DONE
