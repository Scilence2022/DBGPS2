# External analyses: code and summary evidence

Selected files from the working snapshots of the analyses in the manuscript that do
not come from the thermal-aging pipeline. Contents were copied unchanged; paths inside
scripts still point to the original working directories.

| Directory | Contents | Source |
|---|---|---|
| `benchmarks/reconstruction_speed` | matched photolithographic timing harness and tables | `657c8ee:output/benchmarks/2026-10-19` |
| `benchmarks/speed_extension` | format callers, eight-worker search and Bloom matrix | same |
| `benchmarks/resource_control` | serial and three-primary-worker adapters, CPU records | same |
| `simulations/reconstruction_boundary` | exploratory simulations | `657c8ee:output/exploration/2026-10-20` |
| `simulations/boundary_deep_validation` | locked validation, bounds and policy results | same |
| `photolithographic` | adapters, decoding reports, File 1-3 evidence | `9eef8d1c:output/antkowiak2020_reports` |
| `organick` | adapter, validation, optimisation evidence | `4a7c7f2a:experiments/nbt4079/deliverables/2026-10-18/build/DBGPS2_nbt4079_code` |

Short commit ids refer to `657c8ee` and to `refs/codex/snapshots/9eef8d1c65383b026ce6c2e8617545fd8d9e89bd`
and `refs/codex/snapshots/4a7c7f2a833209156dc5eb5d1a3147fec66c46e7` in the local repository.

Not included: raw reads, per-run output directories, decoded archives, compiled files and
third-party source packages. Tables of 2 MB or more are listed in `skipped_large_files.tsv`.
The DNA Fountain adapter, read-merging comparator and decoder live on branch
`dnafountain-adapter`; `../verify_fountain_design.py` scores their outputs.
