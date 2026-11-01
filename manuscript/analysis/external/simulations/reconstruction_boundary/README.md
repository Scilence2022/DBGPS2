# Reconstruction-boundary exploration

Start with `REPORT.md` (Chinese findings) and `PROTOCOL.md` (experimental design and definitions). This is an exploratory mechanism study, not a speed benchmark or a complete-file recovery experiment. No manuscript Word document is modified by these scripts.

## Contents

- `all_simulation_results.csv`: 43,776 parameter points from 1,368 independently counted graphs. Within a condition, three random design pools are the replicates; parameter points are not independent replicates.
- `HT176_exact_boundaries.csv`: 16 parameter points on the original, unmerged paired HT176 reads.
- `runs/`: per-graph settings, generated-input/reference hashes, commands, counters, native checks, and compact count diagnostics. Synthetic reads can be regenerated from the locked settings and seeds.
- `*.truth.u16`, `*.alts.u16`, and associated `*.oracle.json`: lossless reference-path and alternative-branch counts for diagnostic reanalysis. The metadata gives array dimensions. These files are little-endian unsigned 16-bit integers; the counter saturates at 16,383.
- `src/`: frozen production sources plus the diagnostic extension and its original profiler. `bin/boundary_profile` is the measured macOS arm64 binary.
- `reference_evaluation_only/HT176_ref.fa`: known HT176 designs used for scoring and diagnostic ceilings. They do not prune or repair the search.
- `validation_reference/HT176.grid.tsv`: the earlier full-raw-data profile used for 12 overlapping checks.
- `input_revalidation.json`: fresh size/SHA256 checks of both raw inputs and the reference against the previous input provenance.
- `archive_manifest.json`: size and SHA256 of every other archived file; the archive validation sidecar also records the ZIP checksum.

Development inputs, the aborted fixed-terminal-look-ahead diagnostic runs, scratch files, Python caches, and the multi-gigabyte raw FASTQ files are excluded from this delivery. The excluded runs remain isolated in the original workspace and are not in any result table.

## Reanalysis

With Python 3, NumPy, pandas, and Matplotlib installed, run from this folder:

```sh
python3 analyze_results.py
python3 write_report.py
```

These commands regenerate the summary tables, figures, and report from the saved output. They verify the recorded measurement-source/binary hashes; they do not rerun counting or require raw FASTQ files. Figure typography can vary with the installed Matplotlib/font versions. The delivered figures were inspected visually; rerendered figures require inspection again.

## Rerunning measurements

Preserve this archive as the immutable original and use a separate copy for new runs. The measured sources and `execution_manifest.json` must not be silently replaced. To build the diagnostic executable from the saved sources on a compatible system:

```sh
cc -O3 -Isrc -o bin/boundary_profile src/boundary_profile.c src/dbgps_codec.c src/dbgps_redundancy.c src/dbgps_strand.c src/kthread.c -lz -lpthread
```

Run `smoke.py` before measurement. The main simulation is in `run_exploration.py`; the additional, sequential geometry runs are in `geometric_stress.py` and `regular_geometry_control.py`. The runners resume existing case JSON files, so a fresh measurement needs a new output folder without previous run records. Preserve locked designs, record the new executable/source hashes and environment as a separate provenance record, and do not present them as the original execution manifest.

The raw HT176 and historical validation paths in `run_exploration.py` are absolute paths from the measurement machine. Configure both raw FASTQ paths, the evaluation-reference path, and the historical profile path on another machine before that stage. Original raw-file checksums are in `input_revalidation.json`. No public accession for these local files has been verified here. This is a source-and-data record, not a claim of a fully portable one-command workflow.

## Figure legends

1. **Fig_boundary_phase_maps:** correct reconstruction / output-path coverage ceiling (%), at c=1 and r=2.55. Each cell is a mean of three independently generated pools. Rows vary k and internal-fragment retention q; columns vary target retained-base depth. Cut spacing is the regular-partition model parameter, not the observed mean fragment length. The 200* intact control has no q dependence and is reused explicitly. Full results include the omitted 10× depth and q=0.25 conditions.
2. **Fig_boundary_ratio_window:** (a) true-path rule admissibility and unique correct recovery in raw HT176 at c=1; off denotes the native 10^6 ratio ablation, not the published original DBGPS method. (b) terminal-support, ratio-pruning, and search/uniqueness losses at selected r. (c) a synthetic counterexample to the universal optimality of r=2.55: cut spacing 45 nt, target depth 100×, q=0.05, substitution probability 1%, k=21, c=1. (d) matched target-depth q comparison for that synthetic setting. Synthetic curves/bars are means of three pools; all individual results and between-pool standard deviations are in the tables. HT176 is one real library.
3. **Fig_boundary_fragment_geometry:** (a,b) matched regular partitions and Bernoulli cuts on the same three 256-design pools, at target depth 100×, q=1, substitution probability 1%, k=31. Black ticks in (a) mark coverage ceilings; error bars in (a,b) show sample standard deviations across pools. (c) default-parameter recovery versus coverage ceiling for all simulated graphs with nonzero coverage ceiling, colored by mean true-path coverage CV. These are descriptive, correlated conditions, not independent regression samples.

The coverage ceiling here is reference-assisted and specific to output-defining k-mer support, cutoff, seed/terminal conventions, and this implementation. It is not a universal information-theoretic or file-decoding bound. See the report for the full nested decomposition and limitations.
