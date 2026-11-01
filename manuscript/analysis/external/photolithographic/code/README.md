# Antkowiak 2020 reconstruction benchmark

The completed local report and all results are under `input/antkowiak2020/results/`.
The input directory is ignored by Git. Production DBGPS2 sources are unchanged.
No commits or pushes are part of this experiment.

This is an **experimental DBGPS2 adapter**, not a claim that the stock CLI supports
the original synthesis format. `graph_adapter.c` includes DBGPS2.c to reuse its
exact counter, canonical graph lookups, saturation, and coverage-ratio rules.
It enumerates observed k-mer seeds containing the known internal 12-base index,
extends in both directions, and accepts exactly one RS(20,18)/GF(64)-valid 60-base
strand. The 1050-path limit is retained; an incomplete search returns an erasure.
Reference payloads are never used to seed or select a graph path. The C program
reads the optional reference FASTA only after assembly, for graph-support diagnostics.

The numerical index domain and wire format are encoding metadata. They are not
learned per read. The physical format is `42 nt payload | 12 nt scrambled index |
6 nt inner parity`, with index XOR 1071068 and GF(64) polynomial 91.

## Public sources

- Article: https://doi.org/10.1038/s41467-020-19148-3
- Primary code and File 1 design/original ZIP:
  https://github.com/MLI-lab/noisy_dna_data_storage
  at `41422be448b3c1da9041c3a2967d2754d2ad8b1a`.
- RS implementation: https://github.com/reinhardh/dna_rs_coding,
  full checkout at `455e1a5182bf64281f54d434a4ed5ca456c2373e`.
  `build.py` extracts the pre-randomizer-change headers at `8e7e3b5`, strips
  unused Boost dependencies, selects the existing spectral outer decoder,
  and explicitly reproduces the historical libc++ uniform sampler so that
  Linux libstdc++ and macOS libc++ decode the same synthesis randomization.
- Later File 2/3 design deposits:
  https://github.com/fml-ethz/dt4dds-challenges_notebooks
  at `861d3b95b52d1b121e9d85040f20c8ecf0cda2d4`.
- FASTQ file IDs: 22305708, 22305990, and 49503285.
  The original Figshare article IDs are 12128994, 12129090, and 27134316.

`prepare_references.py` is essential for File 2/3: the later deposited designs
are not the synthesized wire layout. Reversing their base order reveals an
unrandomized front index and a 42-base payload. Preserve that public payload,
reserialize the physical scrambled internal index, and regenerate the inner
parity. The conversion reproduces every primary File 1 strand, and File 2/3
then decode into ZIP files with valid internal CRCs and match physical reads.
Keep the original design files and the transcoding manifest for provenance.
Do not score the later design files directly against the physical reads.

## Reproduce

Run commands from the repository root. C/C++ compilers, zlib, Python 3.9+, and
MUSCLE **3.8.31** are required. Matplotlib is needed only for report plots.

```sh
python3 benchmarks/antkowiak2020/fetch_sources.py
ANTK_CHUNK_MB=1 ANTK_DOWNLOAD_WORKERS=64 python3 benchmarks/antkowiak2020/download.py
python3 benchmarks/antkowiak2020/build.py
python3 benchmarks/antkowiak2020/build_muscle.py
python3 benchmarks/antkowiak2020/prepare_references.py

input/antkowiak2020/bin/prepare_reads input/antkowiak2020/raw/I16_S2_R1_001.fastq.gz input/antkowiak2020/processed/file1 > input/antkowiak2020/results/file1_reads.json
input/antkowiak2020/bin/prepare_reads input/antkowiak2020/raw/I18_S3_R1_001.fastq.gz input/antkowiak2020/processed/file2 > input/antkowiak2020/results/file2_reads.json
input/antkowiak2020/bin/prepare_reads input/antkowiak2020/raw/File3.fq.gz input/antkowiak2020/processed/file3 > input/antkowiak2020/results/file3_reads.json

python3 benchmarks/antkowiak2020/run_graphs.py file1
python3 benchmarks/antkowiak2020/run_graphs.py file1 --k 13 15 17 19 21 25 31 --cutoffs 9,19,49,99,199
python3 benchmarks/antkowiak2020/run_graphs.py file1 --k 15 17 19 21 25 31 --cutoffs 1,3,4,5,6
python3 benchmarks/antkowiak2020/run_graphs.py file2 --k 13 15 17 19 21 25 31 --cutoffs 1,3,9,19,49 --parallel 2
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 13 15 17 19 21 25 31 --cutoffs 1,3,9,19,49 --parallel 2

# Follow-up: reuse each k graph across multiple coverage ratios.
cp input/antkowiak2020/results/file3_graph_metrics.json input/antkowiak2020/results/file3_before_ratio_sweep.json
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 15 17 19 21 --cutoffs 1,3,9,19 --ratio 1.2,1.5,1.8,2.1,3.5 --parallel 3
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 13 --cutoffs 9,49,99 --ratio 1.2,1.5,1.8 --parallel 1
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 18 20 --cutoffs 1,3,9 --ratio 1.8,2.1,2.55,3,3.5,5,10 --parallel 1
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 19 --cutoffs 1,3,9 --ratio 3,5,10 --parallel 1
python3 benchmarks/antkowiak2020/run_graphs.py file3 --k 18 --cutoffs 2 --ratio 1.8,2.1,2.55,3.5 --parallel 1
python3 benchmarks/antkowiak2020/file3_ratio_report.py
```

`-c C` requires k-mer count **greater than C**, matching DBGPS2. Ratios use its
integer rounding and saturated counts. Large datasets use the existing two-pass
singleton filter (`ANTK_BLOOM=1` selects a 31-bit bitmap exponent). Removing
count-one keys preserves every query eligible in this grid (all cutoffs >=1).
Early runs used exponent 27; eligibility and counts >=2 are the same.
The follow-up report compares against the preserved original File 3 metric
snapshot `file3_before_ratio_sweep.json`. Ratios in the adapter accept a comma-
separated list. Search output names include each ratio, and run logs include the
entire ratio list, preserving the earlier fixed-ratio experiment logs.

### Candidate-level inner RS correction

The earlier graph output already passes the inner parity check; running an RS
decoder again on that output cannot repair its undetected errors. Revisit every
complete, non-overflow graph path instead. `build_inner_graph.py` generates an
ignored copy of the strict adapter with candidate callbacks, links
`graph_inner_rs.cpp`, and leaves the original adapter and production sources
unchanged. The collector corrects at most one GF(64) symbol, reassigns repaired
indices, and writes decisions before reading an optional evaluation reference.

```sh
python3 benchmarks/antkowiak2020/prepare_inner_grid.py
python3 benchmarks/antkowiak2020/build_inner_graph.py
c++ -O3 -std=c++11 -Ibenchmarks/antkowiak2020 benchmarks/antkowiak2020/test_inner_graph.cpp benchmarks/antkowiak2020/graph_inner_rs.cpp -o input/antkowiak2020/bin/test_inner_graph
input/antkowiak2020/bin/test_inner_graph input/antkowiak2020/sources/File1_design.fasta input/antkowiak2020/results/inner_rs/decision_fixture
python3 benchmarks/antkowiak2020/run_inner_graph.py --parallel 3
python3 benchmarks/antkowiak2020/report_inner_graph.py
```

Use the bundled Python with the installed Matplotlib if adding plots. The inner
report itself only needs the standard library; `--plot` adds a comparison figure.
The completed inner report is `input/antkowiak2020/results/inner_rs/REPORT.zh-CN.md`.
It distinguishes reference fragments appearing at least once among decoded
candidates from correct one-sequence-per-index output decisions. Of 179 original File 3 settings,
119 have fewer complete paths than the attained strict lower bound of 137662
correct strands. One path can map to at most one radius-one RS codeword, so those
settings cannot improve the maximum for a **single configuration**. All remaining
60 settings are rerun; their strict outputs must match the original bytes.
Candidate coverage found using the reference is reported separately and never
presented as a realizable blind selection result. Path votes are not read votes.

### Per-index union across parameter settings

The single-configuration path bound cannot prune settings when measuring a
union: a weak setting may still recover an index missed by stronger settings.
Run all missing original settings, including k=13, without that pruning:

```sh
python3 benchmarks/antkowiak2020/run_inner_graph.py --remaining --k 19 20 18 17 15 13 21 25 31 --parallel 3
python3 benchmarks/antkowiak2020/report_parameter_union.py
python3 benchmarks/antkowiak2020/validate_parameter_union.py
```

`--remaining` resumes completed settings and preserves the earlier run metadata.
The union report validates all 179 settings against the original strict outputs.
It writes `results/inner_rs/parameter_union/REPORT.zh-CN.md`, a per-index table,
a compressed index-by-parameter success matrix, the successful index/parameter
pairs with exact zero/one-symbol path counts, and per-policy success bitmasks.
Parameter IDs refer to `parameters.tsv`; indices are zero-based **corrected
decoded indices**, so repairing an index symbol can reassign a candidate.

A fragment counts once if any tested setting has its exact complete 60-nt
reference after inner RS decoding. Candidate-pool coverage and correct final
one-sequence outputs under each fixed policy are reported separately. Both
cross-parameter success unions use the reference to label successes, and are
not blind candidate selection algorithms or verified file recovery. Missing
candidate counts per block are only an oracle erasure feasibility diagnostic.
The earlier single-configuration report is retained as the historical
60-configuration result; the new union report includes all 179 configurations.

### DBGPS-analyzer diagnosis of block 12

`run_block12_analyzer.py` invokes the user's checkout at
`/Users/song/.codex/worktrees/dabb/DBGPS-analyzer` without changing that repository.
Build `DBGPS-analyzer.c` plus `kthread.c` and `vendor/bfc/bbf.c` with
`cc -g -Wall -Wextra -O3 -lz -lpthread`, and build `DBGPS-links.c` plus
`kthread.c` with the same flags. Place the two binaries in
`input/antkowiak2020/bin/analyzer_tools/` before running:

```sh
python3 benchmarks/antkowiak2020/run_block12_analyzer.py --parallel 3
PYTHONPATH=input/antkowiak2020/plotlibs /Users/song/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 benchmarks/antkowiak2020/report_block12_analyzer.py
```

The report needs NumPy and Matplotlib; these were already installed in the
ignored local plotting directory. It compares all 12 blocks at the nine
previously tested k values, with block 12 as the focus. Target strands retain
all 60 nt with zero primer trimming. All retained 92,337,721 processed reads
and all 195,619,515 original reads are used; no subsampling is involved.

There are 13 analyzer runs: processed reads at all nine k values; original
reads at k=17/19/20 with actual lengths preserved; and original reads at k=19
capped at 60 nt. A full read index for processed k=19 measures global noise;
the other runs use the exact-reference Bloom backend, which preserves target
counts but does not measure off-reference noise. DBGPS-links analyzes the
complete design and each block separately at all nine k values.

The standalone coverage output deduplicates canonical keys within each strand.
`report_block12_analyzer.py` expands repeated keys back to physical positions
before analyzing profiles. It checks every processed minimum against the old
filtered counter, validates the CLI's Sm/unique-target Kd totals, and predicts
the exact unmodified reference path under all 179 directed integer-rounded
rules. After accounting for original overflow statuses, those predictions
must equal the observed zero-error-candidate flags at every index.

Results, raw analyzer outputs, commands, provenance and SHA-256 checks are in
`input/antkowiak2020/results/analyzer_block12/`. The report distinguishes
reference-path coverage from actual reconstruction, pooled k-mer counts from
assigned read depth, and the current input's RS recovery from potential gains
under a different preprocessing scheme. No new reconstruction from raw reads
is claimed by this diagnosis.

### Align reads to the design, then decode blocks independently

This is explicitly **reference-assisted partitioning**: minimap2 uses the full
known design payload to assign raw reads to blocks. It is not a reference-free
end-to-end recovery benchmark. The splitter only routes full original reads;
it does not replace bases, fill gaps, or create consensus from the design.
Reverse-strand SAM primary sequences/qualities are restored to their original
orientation before FASTQ export. Unassigned and cross-block ambiguous reads
are retained. Both broad and strict files are exported for every block.

```sh
c++ -O3 -std=c++11 -Wall -Wextra -Werror benchmarks/antkowiak2020/split_aligned_blocks.cpp -lz -o input/antkowiak2020/bin/split_aligned_blocks
python3 benchmarks/antkowiak2020/build_block_graph.py
python3 benchmarks/antkowiak2020/split_reads_by_alignment.py
python3 benchmarks/antkowiak2020/run_block_decode.py --phase controls --parallel 2
python3 benchmarks/antkowiak2020/run_block_decode.py --phase assigned-control --parallel 2
python3 benchmarks/antkowiak2020/run_block_decode.py --phase blocks --parallel 3
PYTHONPATH=input/antkowiak2020/plotlibs /Users/song/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 benchmarks/antkowiak2020/report_block_decode.py
```

All 195,619,515 raw reads are aligned; SAM streams directly to the C++ router,
so a large full SAM is not stored. `split_provenance.json` records commands,
thresholds and counts. The local first-100,000-read pilot verified the original
sequence and quality for every exported record and complete disjoint routing;
its measured evidence is retained in `pilot_split/sequence_validation.json`.
The report reads that evidence file when available rather than manufacturing
a pilot pass on a fresh run.

`graph_block_rs` is a separate generated binary whose only search-program
change is a configurable read-counting length cap. At limit 60 all eleven
decision outputs must exactly reproduce the original k19/c1/r3.5 result.
New raw-read runs use limit 200 (the raw reads are no longer truncated), the
same 1050-live-path limit, full global address domain, and unchanged inner RS
collector. Do not renumber a block's indices to 0..16382: its encoded addresses
remain global. The reference is still loaded only after output decisions by
the graph program; the earlier mapping step is separately reference-assisted.

The matched grid has ten settings: k15/c1/r1.8,2.55; k17/18/19/c1/r2.55,3.5;
and k20/c1/r3.5,5.0. Every broad block and the exact same assigned reads
rejoined into one pool use that grid. The rejoined-pool control
isolates the effect of separating blocks from alignment-based read filtering.
Block 12 broad additionally uses all 179 original settings; its strict and
ambiguity-inclusive variants use the matched ten. There are 319 new scored
configurations plus the byte-identical length-60 control.

Unfiltered all-raw mixed controls are optional with
`--phase controls --include-unfiltered-raw`. The present experiment retained
five completed auxiliary settings, then stopped the remaining unfiltered-raw
work after k18/c1/r2.55 took 700.778 search seconds and overflowed
196322/196596 indices. Incomplete settings are not counted as zero and are
excluded from the matched comparison. The key rejoined-versus-separated
assigned-read comparison still completes all ten settings.

`audit_fastq_multiset.cpp` uses macOS CommonCrypto SHA-256 record digests to
check the complete source and exported name/sequence/quality multisets. The
full 195,619,515-read check and the separate-versus-rejoined assigned-pool
check are recorded in `full_fastq_audit.json` and `assigned_pool_audit.json`.

The report distinguishes candidate unions from actual unique decisions and
computes six-lane outer-RS budgets. If actual outputs for all twelve blocks
meet the bound, it attempts a real outer decode and validates the resulting
ZIP CRC and bytes. Large FASTQs and all outputs remain ignored under
`input/antkowiak2020/results/block_decode/`.

`build_muscle.py` verifies and builds the original MUSCLE 3.8.31 into
`input/antkowiak2020/bin/muscle3`. On Apple Silicon it removes the old build script's
x86 `-msse2 -mfpmath=sse` flags and compile its listed translation units with
`clang++ -O3 -DNDEBUG=1 -D_FILE_OFFSET_BITS=64`. No alignment parameters were changed.

```sh
python3 benchmarks/antkowiak2020/msa_baseline.py input/antkowiak2020/processed/file1.fa.gz input/antkowiak2020/sources/File1_design.fasta input/antkowiak2020/bin/muscle3 input/antkowiak2020/results/file1_msa
python3 benchmarks/antkowiak2020/msa_baseline.py input/antkowiak2020/processed/file2.fa.gz input/antkowiak2020/sources/File2_synthesized.fasta input/antkowiak2020/bin/muscle3 input/antkowiak2020/results/file2_msa --workers 8
python3 benchmarks/antkowiak2020/msa_baseline.py input/antkowiak2020/processed/file3.fa.gz input/antkowiak2020/sources/File3_synthesized.fasta input/antkowiak2020/bin/muscle3 input/antkowiak2020/results/file3_msa --workers 16
python3 benchmarks/antkowiak2020/msa_baseline.py input/antkowiak2020/raw/I16_S2_R1_001.fastq.gz input/antkowiak2020/sources/File1_design.fasta input/antkowiak2020/bin/muscle3 input/antkowiak2020/results/file1_notebook --workers 8 --prefix-length 14 --min-cluster 4 --gap-weight 0.5 --raw-fastq
python3 benchmarks/antkowiak2020/report.py
```

The matched baseline follows the paper's prefix16/minimum5/gap0.4 settings on
the same length-filtered and capped reads as the graph. The additional File 1
run follows the released notebook's prefix14/minimum4/gap0.5 and aligns before
capping. Both use five random samples of 15 for larger clusters, fixed seed
20261017, and the full retained datasets. No full LSH reproduction is claimed.
`inner_select` adds a common inner-code/voting comparison, including a separately
reported one-symbol-correcting variant; ties become erasures.

`rs_driver_fast` uses the authors' decoder and an equivalent GF(2^14) log table
for products, inverses, and powers. `test_fast_gf.cpp` compares it with the
original bitwise field; original and accelerated File 1 outputs are also
checked byte-for-byte against the author's ZIP. A decoded reference ZIP is a
positive control; a `fileN_best.zip` is reconstructed from sequencing candidates.
The report distinguishes these artifacts and verifies their bytes and ZIP CRCs.

Timing logs are collected under concurrent workloads and different worker
counts. Use them as local resource records, not as controlled speed benchmarks.
The best parameters were selected using these datasets' truth, without an
independent holdout; the parameter envelope is exploratory.

For local plotting dependencies, use `uv pip install --python python3 --target
input/antkowiak2020/plotlibs matplotlib`, then run `report.py` with that same
Python version. The current report uses Codex's bundled Python runtime.
