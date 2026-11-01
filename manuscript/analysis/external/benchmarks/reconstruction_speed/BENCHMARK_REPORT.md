# Reconstruction speed comparison

Mean ± sample standard deviation of three serial computational repeats. Inputs are the same processed gzipped FASTA for every method. Exact recovery and wrong counts refer to unique strict strand outputs, before outer decoding.

| Input | Method | Time (s) | Exact recovery (%) | Wrong outputs | Speedup vs MSA |
|---|---|---:|---:|---:|---:|
| file1 | Clustering + MSA + vote | 267.729 ± 4.857 | 65.7877 | 4 | 1.000× |
| file1 | DBGPS rule | 14.776 ± 0.339 | 93.8351 | 0 | 18.119× |
| file1 | DBGPS2 | 15.894 ± 0.744 | 98.4252 | 0 | 16.845× |
| file2 | Clustering + MSA + vote | 483.451 ± 4.444 | 22.9181 | 13 | 1.000× |
| file2 | DBGPS rule | 26.923 ± 0.120 | 79.5601 | 11 | 17.957× |
| file2 | DBGPS2 | 34.020 ± 1.282 | 88.4494 | 3 | 14.211× |

The DBGPS comparator is its absolute-cutoff rule emulated in the same format adapter, with a pilot-selected cutoff. DBGPS2 uses the existing c=1, r=2.55 settings. The prefix16/MUSCLE 3.8.31 baseline is not the original full LSH pipeline. These are implementation-specific operating-point comparisons; recovery levels differ. The baseline invokes MUSCLE for each cluster/subset, including launch and communication costs.

See README.md for the complete protocol, exclusions, configured worker counts and limitations. See all_repeats.csv, stage_timings.csv, raw run logs and output SHA256 digests for the evidence.
