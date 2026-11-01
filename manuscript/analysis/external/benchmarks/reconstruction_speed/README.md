# Controlled reconstruction timing

This benchmark compares the adapted DBGPS absolute-coverage rule, adapted
DBGPS2 local-ratio search, and a conventional prefix-clustering, MUSCLE 3.8.31,
weighted-majority-voting pipeline on the same preprocessed photolithographic
inputs. It is an implementation benchmark, not a ranking of all clustering or
trace-reconstruction methods.

## Inputs and experimental scope

- File 1: 14,299,688 retained reads, 16,383 designed 60-nt codewords.
- File 2: 33,404,713 retained reads, 49,149 designed 60-nt codewords.
- Inputs are the existing `file1.fa.gz` and `file2.fa.gz` from the Antkowiak
  reanalysis. Every method reads the same gzipped FASTA for a dataset.
- The common filtering and truncation that produced these FASTAs is excluded
  from every timing. Decompression and input parsing are included.
- This experiment tests highly erroneous photolithographic libraries. It does
  not measure the speed of MSA reconstruction on the thermal-aging samples.
- Reconstruction never loads a design sequence. Post-run evaluation compares
  unique strict outputs with the corresponding physical-format designs.

## Methods

The graph adapter uses the exact DBGPS2 counter and internal-index seed
enumeration, left/right traversal and inner RS validity test previously used
in this reanalysis. Both graph modes have the same 1,050-live-path limit,
single-threaded search and eight hash-insertion workers for counting. The
counter additionally pipelines input reading and k-mer extraction, so eight
is the configured worker count rather than a hard cap on all active threads.
Both Bloom filters are
disabled. Only the local-ratio rule and absolute cutoff differ. `r=1000000`
emulates the DBGPS absolute-cutoff rule; it is not the originally published
2022 executable. The seed count threshold in this adapter excludes singleton
seeds in both modes. File 1 uses k=19, File 2 k=21.

DBGPS2 uses c=1 and r=2.55, fixed from the existing reanalysis. The DBGPS
cutoff was selected retrospectively by exact recovery in the pilot grid
c={1,3,5,8,10,15,20,30}, using the same k. The selected cutoffs are 5 and 3.
These are exploratory, dataset-tuned comparisons. Pilot costs are saved but
excluded from the repeated timings; this is not time to discover parameters
for an unseen library.

The conventional baseline clusters by the exact first 16 nt, retains clusters
with at least five reads, and runs MUSCLE 3.8.31 followed by per-column voting.
Clusters of at most 15 reads are aligned in full; larger clusters produce
five subsets of 15 reads, sampled using the fixed seed 20261017 and prefix.
Gaps receive vote weight 0.4. Only 60-nt ACGT candidates proceed to the
existing inner-code selector; strict outputs require a zero-syndrome codeword
and an untied vote leader per index. The selector also writes its corrected
output, but only the strict output is evaluated here. Eight independent
single-threaded alignment workers are used. This is the previously reported
prefix baseline, not a reproduction of the original full LSH pipeline and
not the alternative notebook-style prefix14 baseline.

## Timing protocol

All 18 measured processes run serially on the same Apple M3 Max workstation
(16 CPU cores, 128 GiB physical RAM), with three repeats per method and input.
Method order rotates between repeats and datasets. Background user
applications remain present; process snapshots before each run document the
load. The benchmark launches no concurrent analysis jobs. No filesystem cache
is forcibly dropped. Results describe warm/mixed-cache workstation use.
Thread environment variables cap numerical-library threading at one.

End-to-end elapsed time includes process startup, input reading, clustering or
graph construction, candidate generation, inner-code decisions, output writing
and process teardown. Post-run truth scoring, pilot selection, downloading,
shared preprocessing, and outer file decoding are excluded. `/usr/bin/time -l`
is retained in each stderr log. Graph counters use eight insertion workers
plus reading/extraction pipeline stages; graph search is serial. The baseline
uses eight alignment workers plus coordinator overhead. These are matching
configured worker counts, not identical total thread counts, CPU affinity or
realized CPU utilization at every stage. The MSA/voting phase is parallel and its wall
time is reported jointly. Per-worker elapsed sums are diagnostic and must not
be mistaken for end-to-end wall time or CPU time.

Report mean and sample standard deviation over three computational repeats,
alongside exact recovery and erroneous outputs. Speedup is the ratio of mean
MSA-pipeline time to mean graph-pipeline time at the reported operating points.
Recovery rates differ; these are not comparisons at matched recovery or file
decoding endpoints. They do not establish universal algorithmic complexity
advantages, gains attributable solely to local ratio pruning, or speedups over
every clustering/MSA implementation.

## Reproduction

Use the bundled Python runtime in the recorded commands. Build the two C/C++
tools from the source snapshot with Apple clang 17.0 and `-O3`. The MUSCLE
binary reports v3.8.31; its source archive is retained in the original input
directory and is available from the author's legacy release site:
https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_src.tar.gz

`requirements.txt` records the Python versions used for plotting and document
generation. Reconstruction itself uses the standard library and the retained
executables. Figures and tables can be regenerated from `results.json` without
the original sequencing reads. QA page renders are not included in this package.

1. Use `build_tools.py` for the counter, adapter, selector and optional read
   preparation tool. The retained MUSCLE binary is macOS arm64; use
   `build_muscle.py` to rebuild it from the included hash-checked source archive
   for another supported architecture. Compilation is outside timed runs.
2. If recreating the processed input, `prepare_reads original.fastq.gz prefix`
   retains 55–70-nt reads and caps sequences at 60 nt. Original FASTQ source
   file identifiers and checksums are recorded in the manuscript and input
   provenance. Evaluation references are included under
   `reference_evaluation_only/` and must never be passed to reconstruction.
3. Regenerate the pilot graph outputs with the saved grid and design-free
   reconstruction argument `-`.
4. Run `run_benchmark.py` with the input root adjusted if needed. It resumes
   completed runs from their result records.
5. Run `validate_results.py` and `summarize_and_plot.py` after all 18 runs
   complete. Historical compatibility validation requires the retained prior
   baseline outputs; its path can be adjusted when reproducing elsewhere.

For a new input location, set `DATA` in `run_benchmark.py` to a directory with
`processed/file1.fa.gz`, `processed/file2.fa.gz` and `sources/` containing the
two included evaluation references. Pilot reconstruction commands are:

```text
bin/graph_adapter processed/file1.fa.gz 16383 19 pilot/file1 - 1000000 1,3,5,8,10,15,20,30
bin/graph_adapter processed/file2.fa.gz 49149 21 pilot/file2 - 1000000 1,3,5,8,10,15,20,30
```

`source_sha256.json` pins the core, adapter and executable snapshot. The adapter
diff from its source is diagnostic timing only: the seed-enumeration log adds
elapsed time from the start. The reconstruction logic is unchanged. Results
contain exact commands, raw elapsed times, stage timing, exact output counts,
incorrect counts, and SHA256 digests of strict outputs. Summary generation
requires identical strict-output hashes across the three repeats.
