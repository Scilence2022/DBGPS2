#!/bin/bash
# Log system CPU usage and load every 10 s so that each timing run can be checked for interference.
while [ ! -f /Users/song/dbgps2_ht/results/TIMING_CLEAN_DONE ]; do
  echo "$(date +%H:%M:%S) $(uptime | sed 's/.*load averages: //') $(top -l 1 -n 0 | grep 'CPU usage' | sed 's/CPU usage: //')"
  sleep 10
done
