#!/bin/bash
# usage: subsample.sh in.fq.gz out_prefix every
f=$1; o=$2; n=$3
gzip -dc "$f" | awk -v n=$n -v o="$o.fq" 'BEGIN{r=0} {l=NR%4; if(l==1){r++; keep=(r%n==0)} if(keep) print > o} END{print r > (o ".count")}'
