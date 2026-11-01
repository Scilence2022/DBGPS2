#!/bin/bash
# Sequencing-depth series: nested systematic subsamples (read pair j kept if (j-1) mod 100 < p).
# 50% -> 20% -> 10% files are materialised; 1/2/5% are streamed from the 10% file.
cd /Users/song/dbgps2_ht
W=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19
O=results/depth; mkdir -p $O fifo fq50 fq20 fq10
for f in fq/*.fq; do b=$(basename $f)
  [ -s fq50/$b ] || (awk 'NR%4==1{r++; keep=((r-1)%100<50)} keep' $f > fq50/$b) &
done; wait
for f in fq50/*.fq; do b=$(basename $f)
  (awk 'NR%4==1{r++; keep=((r-1)%50<20)} keep' $f > fq20/$b && \
   awk 'NR%4==1{r++; keep=((r-1)%20<10)} keep' fq20/$b > fq10/$b) &
done; wait
run() { s=$1; K=$2; k=$3; tag=${s}_k${k}_p$K
  case $K in
    50|20|10) p1=fq$K/${s}_1.fq; p2=fq$K/${s}_2.fq;;
    *) p1=fifo/${tag}_1.fq; p2=fifo/${tag}_2.fq; rm -f $p1 $p2; mkfifo $p1 $p2
       for r in 1 2; do
         (awk -v K=$K 'NR%4==1{r++; keep=((r-1)%10<K)} keep' fq10/${s}_$r.fq > fifo/${tag}_$r.fq) &
       done;;
  esac
  ./dbgps2_profile_v2 -N -k $k -x ref.fa -M $W/examples/legacy-2base.map \
    -c 1,2,3,4,5,6,8,10,12,15,20,30,40 -r 0,2.55,3 -t 4 -o $O/$tag $p1 $p2 > $O/$tag.stdout 2> $O/$tag.stderr
  rm -f fifo/${tag}_1.fq fifo/${tag}_2.fq; }
for k in 31 21; do for K in 1 2 5 10 20 50; do for s in HT_con HT1 HT2 HT3 HT4 HT5; do
  while [ $(pgrep -f "dbgps2_profile_v2 -N" | wc -l) -ge 4 ]; do sleep 3; done
  run $s $K $k &
  sleep 1
done; done; done
wait
echo done > $O/ALL_DONE
