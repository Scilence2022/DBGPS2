#!/bin/bash
# HT176 (176 days at 70 C): full-depth decoding grid, k-mer size sweep and sequencing-depth series.
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
MAP=$W/examples/legacy-2base.map
CS=1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100
s=HT176
mkdir -p results/full results/ksweep results/depth fifo fq50 fq20 fq10
# 1) full depth, k = 31
./dbgps2_profile_v2 -x ref.fa -M $MAP -c $CS -r 0,1.5,2,2.55,3,4,6 -t 16 -T trace_ids.txt \
  -o results/full/$s fq/${s}_1.fq fq/${s}_2.fq > results/full/$s.stdout 2> results/full/$s.stderr
# 2) k-mer size sweep, two runs at a time
for k in 21 25 19 23 27 17; do
  while [ $(pgrep -f "dbgps2_profile_v2 -k" | wc -l) -ge 2 ]; do sleep 5; done
  ./dbgps2_profile_v2 -k $k -x ref.fa -M $MAP -c $CS -r 0,2,2.55,3,4 -t 8 \
    -o results/ksweep/${s}_k$k fq/${s}_1.fq fq/${s}_2.fq > results/ksweep/${s}_k$k.stdout 2> results/ksweep/${s}_k$k.stderr &
  sleep 2
done; wait
# 3) depth series: nested systematic subsamples (read pair j kept if (j-1) mod 100 < p)
for r in 1 2; do
  [ -s fq50/${s}_$r.fq ] || awk 'NR%4==1{n++; keep=((n-1)%100<50)} keep' fq/${s}_$r.fq > fq50/${s}_$r.fq
  [ -s fq20/${s}_$r.fq ] || awk 'NR%4==1{n++; keep=((n-1)%50<20)} keep' fq50/${s}_$r.fq > fq20/${s}_$r.fq
  [ -s fq10/${s}_$r.fq ] || awk 'NR%4==1{n++; keep=((n-1)%20<10)} keep' fq20/${s}_$r.fq > fq10/${s}_$r.fq
done
run() { K=$1; k=$2; tag=${s}_k${k}_p$K
  case $K in
    50|20|10) p1=fq$K/${s}_1.fq; p2=fq$K/${s}_2.fq;;
    *) p1=fifo/${tag}_1.fq; p2=fifo/${tag}_2.fq; rm -f $p1 $p2; mkfifo $p1 $p2
       for r in 1 2; do (awk -v K=$K 'NR%4==1{n++; keep=((n-1)%10<K)} keep' fq10/${s}_$r.fq > fifo/${tag}_$r.fq) & done;;
  esac
  ./dbgps2_profile_v2 -N -k $k -x ref.fa -M $MAP -c 1,2,3,4,5,6,8,10,12,15,20,30,40 -r 0,2.55,3 -t 4 \
    -o results/depth/$tag $p1 $p2 > results/depth/$tag.stdout 2> results/depth/$tag.stderr
  rm -f fifo/${tag}_1.fq fifo/${tag}_2.fq; }
for k in 31 21; do for K in 1 2 5 10 20 50; do
  while [ $(pgrep -f "dbgps2_profile_v2 -N" | wc -l) -ge 3 ]; do sleep 3; done
  run $K $k &
  sleep 1
done; done; wait
echo done > results/HT176_DONE
