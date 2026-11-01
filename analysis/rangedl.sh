#!/bin/zsh
# Parallel range download from ENA (server throttles per connection).
# Usage: rangedl.sh <url> <outfile> <size> <nchunks>
set -e
URL=$1
OUT=$2
SIZE=$3
N=$4
CHUNK=$(( (SIZE + N - 1) / N ))
mkdir -p ${OUT}.parts
seq 0 $((N-1)) | xargs -P $N -I{} zsh -c "
  i={}
  START=\$(( i * $CHUNK ))
  END=\$(( START + $CHUNK - 1 ))
  if [ \$END -ge $SIZE ]; then END=\$(( $SIZE - 1 )); fi
  if [ \$START -ge $SIZE ]; then exit 0; fi
  PART=${OUT}.parts/\$(printf %04d \$i)
  for attempt in 1 2 3 4 5; do
    curl -sS --max-time 1800 --retry 3 -r \${START}-\${END} -o \$PART '$URL' && exit 0
    sleep 5
  done
  exit 1
"
cat ${OUT}.parts/* > $OUT
rm -rf ${OUT}.parts
ls -l $OUT
