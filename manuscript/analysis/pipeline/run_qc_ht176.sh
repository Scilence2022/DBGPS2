#!/bin/bash
# Read-level QC of HT176 (176 days at 70 C): same protocol as run_qc.sh (every 25th pair, first 1M reads per mate).
cd /Users/song/dbgps2_ht
A=/Users/song/Github-Repos/DBGPS2/.claude/worktrees/dbgps2-thermal-aging-manuscript-b04b19/manuscript/analysis
mkdir -p results/qc
for r in 1 2; do
  head -n 100000000 fq/HT176_$r.fq | awk 'NR%4==1{n++; keep=(n%25==0)} keep' > results/qc/HT176_$r.fq
  minimap2 -x sr -c --cs=long -t 4 --secondary=no ref.mmi results/qc/HT176_$r.fq > results/qc/HT176_$r.paf 2>/dev/null
  python3 $A/qc_align.py HT176 results/qc/HT176_$r.fq results/qc/HT176_$r.paf results/qc/HT176_$r > results/qc/HT176_$r.summary
  rm results/qc/HT176_$r.fq
done
python3 $A/qc_fragments.py results/qc HT176 > results/qc/HT176.fragments.jsonl
echo done > results/qc/HT176_DONE
