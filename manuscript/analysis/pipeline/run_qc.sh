#!/bin/bash
cd /Users/song/dbgps2_ht
A=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19/manuscript/analysis
[ -f ref.mmi ] || minimap2 -x sr -d ref.mmi ref.fa 2>/dev/null
for s in HT_con HT1 HT2 HT3 HT4 HT5; do for r in 1 2; do
  # first 1M reads of the 4% systematic subsample
  head -n 4000000 sub/${s}_$r.fq > results/qc/${s}_$r.fq
  minimap2 -x sr -c --cs=long -t 4 --secondary=no ref.mmi results/qc/${s}_$r.fq > results/qc/${s}_$r.paf 2>/dev/null
  python3 $A/qc_align.py $s results/qc/${s}_$r.fq results/qc/${s}_$r.paf results/qc/${s}_$r > results/qc/${s}_$r.summary
  rm results/qc/${s}_$r.fq
done; done
python3 $A/qc_fragments.py results/qc HT_con HT1 HT2 HT3 HT4 HT5 > results/qc/fragments.jsonl
echo done > results/qc/ALL_DONE
