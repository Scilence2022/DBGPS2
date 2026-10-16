#!/bin/bash
# DBGPS2 executable benchmarks on a quiet machine: exact vs Bloom options and DBGPS emulation
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/bloom; mkdir -p $O
for s in HT_con HT3 HT5; do
  f1=fq/${s}_1.fq; f2=fq/${s}_2.fq
  for mode in exact F34 B12 F34B12 dbgps; do
    case $mode in
      exact) opt="-r 2.55";; F34) opt="-r 2.55 -F 34";; B12) opt="-r 2.55 -B 12";;
      F34B12) opt="-r 2.55 -F 34 -B 12";; dbgps) opt="-r 1000000 -c 10";;
    esac
    /usr/bin/time -l $W/DBGPS2 -M $W/examples/legacy-2base.map -t 16 -c 1 $opt -o $O/${s}_$mode $f1 $f2 \
      > $O/${s}_$mode.stdout 2> $O/${s}_$mode.time
  done
done
echo done > $O/ALL_DONE
