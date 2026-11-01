# DBGPS2 decoding experiment for Organick et al nbt4079

Prepared 2026-10-18. This is a strand reconstruction experiment, not a completed
outer Reed-Solomon decoder. All reported recovery is before outer RS correction.

## Quick start

Extract into a directory without whitespace; the historical C command protocol separates paths by whitespace.

Use Python 3.10 or newer, a C compiler providing cc, zlib development headers,
and pthreads. The original run used macOS arm64 and Apple clang 17.0.0.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_suite.py --stage check
python run_suite.py --stage smoke
python validate_package.py
python run_suite.py --stage all
```

The smoke command needs no network or large input files. It checks the included
20-strand perfect fixture with two configurations and the substituted-read
fixture with c1. It does not constitute a rerun of the full sequencing dataset.
The added validate_package.py compiles the original-validator, inner-validator and
path-bound drivers, checks both directions at six k values on perfect fixtures,
and checks positive and systematic-gap bound controls. Results are in
ADAPTER_VALIDATION.json. These are packaging checks performed on 18 October 2026,
separate from the archived 17 October experiment controls.
The all command downloads pinned author data and executes all completed phases.
Run individual stages in the order printed by `python run_suite.py --help` to
resume. The original count_stream utility is retained, but the stage wrapper
also creates references.fa, an input needed by the path-bound diagnostic.

## Contents

- Root C sources and headers: frozen source snapshot with LICENSE.
- experiments/nbt4079/scripts: nine original experiment scripts and a provenance helper.
- experiments/nbt4079/optimization/scripts: ten optimization scripts and a provenance helper.
- historical_scripts: byte-for-byte original experiment scripts before packaging.
- fixtures: reference-derived control sequences only; never real-data scoring inputs.
- evidence: original small metrics, all 480 screened configuration scores, nine full
  confirmation records, 19 residual refinement records, ensemble rule grids,
  logs, controls, LFS pointers, and source provenance. These are historical records.
- SOURCE_SNAPSHOT.json and SHA256SUMS: version and integrity records.

The two runnable scripts run_experiment.py and optimize.py differ from their
historical copies only in provenance collection. A standalone ZIP has no Git
history, so the runnable copies read SOURCE_SNAPSHOT.json and compare source
hashes rather than querying git. No decoding logic has been changed.

## Data and disk use

The code archive omits raw reads, all reference designs, large intermediate TSVs,
and prediction FASTAs. It retains the four tiny control fixtures. The two large
author downloads total 1,381,394,524 compressed bytes and are fetched and SHA256
verified by download_data.py. Additional FASTAs, seeds and intermediate TSVs
require several GB; allow at least 15 GB free as a working allowance, not a
measured peak. Prior graph construction used four threads per process, and the
search runner may build two graphs concurrently. No measured peak RAM is available.

Historical evidence contains absolute paths from the original machine. Do not
copy that evidence into the live results directories to resume execution: the
corresponding large candidate files are intentionally absent. A fresh staged run
generates metadata using the extraction path. Evidence is not changed by the runner.

## Analysis details

The raw input contains 15,719,903 reads; the target is 607,150 designed strands.
All full-read optimization graphs use every read. The 1/128 hash sample is a seed
query sample, never a read sample. The final ensemble groups candidates by the
first 13 canonical forward-body bases. Earlier optimization JSON scores use
14-base conflict groups; final_summary.tsv is the definitive 13-base final table.
The baseline summary retains its original 14-base convention. Neither encoded
prefix is a fully decoded block/column address.

References are used for global format calibration, parameter and rule selection,
and exact scoring. They do not supply individual real-data seeds, repair a
candidate, cluster reads, or choose a correct candidate separately for each target.
These are dataset-tuned results, not held-out estimates. Zero observed errors on
the tuning dataset is not an error-free guarantee on new data.

Best exact result: 596,006 correct / 607,150 targets = 98.16454%, plus 166 wrong
outputs. Native-validator tuning with bidirectional search: 595,078 correct.
The inner-code validator and ensemble are extra format adaptation.
The 596,189 path-existence count uses true designs and is diagnostic only.

The optional no-ratio baseline diagnostic was stopped. Its partial files are
excluded from all scores and from this archive. The default runner does not run it.
The authors' BMA decoder was not reproduced. Final digital files and paper-style
cross-strand RS correction have not been validated.

Paper: https://doi.org/10.1038/nbt.4079
Author data: https://github.com/uwmisl/data-nbt17
Author data revision: 1077c74a0e002be83ab832e69bfd560e94d4d202
DBGPS2 source revision: f7b1373ec68bc6518af91b7427eda53699a3c73d
