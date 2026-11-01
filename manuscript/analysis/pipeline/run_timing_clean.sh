#!/bin/bash
# All timing measurements, one at a time, 3 threads, on an idle machine (no other applications; no sleep).
cd /Users/song/dbgps2_ht
./run_bloom2.sh
./run_speed2.sh
./run_speed3.sh
./run_bloom_sweep.sh
echo done > results/TIMING_CLEAN_DONE
