# Reconstruction speed across DNA storage formats and Bloom settings

This package extends the separate photolithographic File 1/2 benchmark. The current matrix uses DNA Fountain, Organick File 16, and HT176. HT5 measurements were excluded when the user replaced that sample with HT176; they do not contribute to the matrix figures or tables. No manuscript files were committed.

## Input identity

- **DNA Fountain:** ENA ERR1816980, 16,692,317 raw read pairs, 72,000 synthesized droplets. NGmerge v0.3 produced 15,589,735 merged molecules; 10,477,842 had the expected 152-nt length and were retained. Reconstruction reports RS(38,36)-zero-syndrome droplets with unique 16-nt seeds. The original complete synthesis list is unavailable here, so this is not an exact-design recovery endpoint. Untimed outer decoding checks archive MD5 `8651e90d3a013178b816b63fdbb94b9b`. Independent re-encoding checks all accepted payloads against the verified archive without inferring membership in the original synthesis seed list.
- **Organick File 16:** Nature Biotechnology (2018), DOI https://doi.org/10.1038/nbt.4079, author data https://github.com/uwmisl/data-nbt17 at revision `1077c74a0e002be83ab832e69bfd560e94d4d202`. The `id20` filenames correspond to File 16: 607,150 designed strands and 15,719,903 raw FASTQ records. The full untrimmed FASTA conversion is used. Exact evaluation compares all 110 encoded bases, with all designs as denominator. No outer RS decoding or digital-file recovery is asserted.
- **HT176:** the local 176-day, 70°C thermal-aging sample of the 210,000-strand Song-format pool. Raw files are `HT176_FDDP220002706-1a_1.fq.gz` and `_2.fq.gz`; sizes, SHA256 and MD5 are in `input_provenance.json`. The sample contains 60,642,860 raw pairs; NGmerge retained 49,620,385 merged molecules. Exact evaluation compares the complete 164-nt index/payload/CRC region over all 210,000 designs. Song et al., Nature Communications 13, 5361 (2022), https://doi.org/10.1038/s41467-022-33046-w, identifies the encoding design and pool format. A public accession or publication of this specific HT176 sequencing sample has not been independently verified; it is not labelled as the paper's public 70-day data. The supplied original article and supplement contain no HT176/176-day designation.

References are not used in reconstruction. They are used in explicitly separate format checks, smoke tests, parameter scoring and output evaluation. Raw reads and large intermediate FASTA files remain local and are excluded from the portable package. Common preprocessing is reference-free. Failed merges are excluded equally from every method; the comparison is conditional on successful merging. It does not replace the full-raw-read thermal experiments elsewhere in the manuscript.

## Nine configurations per dataset

Each of these points has three complete computational timing repeats, for 81 matrix entries:

1. Exact-prefix clustering, MUSCLE 3.8.31, weighted majority voting and format validation.
2. The DBGPS absolute-cutoff rule emulated with r=1,000,000 and no Bloom filters. This measures the rule in the same saved core, not the original published executable.
3. DBGPS2 without Bloom filters at the primary coverage cutoff.
4. DBGPS2 without Bloom filters at an alternative coverage cutoff.
5. DBGPS2 with query filter B8.
6. DBGPS2 with query filter B16.
7. DBGPS2 with two-pass counting filter F32.
8. DBGPS2 with two-pass counting filter F34.
9. DBGPS2 with query filter B12 and two-pass counting filter F34 together.

Primary `(k,c,r)` settings are Fountain `(21,0,6)`, Organick `(29,0,1.9)` and HT176 `(31,1,2.55)`. All five Bloom points retain these primary search settings. Alternative no-Bloom c is 1 for Fountain/Organick and 5 for HT176. The DBGPS rule uses c=1,3,5, respectively. Its cutoffs were selected from a fixed hashed 1/128 seed-sample pilot; the grid was {1,3,5,10,20}, extended with {25,40} for HT176. Pilot selection maximized exact outputs for Organick/HT176 and code-valid count for Fountain. This retrospective tuning is not an independent held-out test and is excluded from timing.

B denotes a minimum query bits-per-key budget; actual power-of-two allocations are recorded in bytes. F is log2 total first-pass bitmap bits: F32 requests 512 MiB, F34 2 GiB. Query filtering preserves exact counts and exact-table fallback. Two-pass counting discards singletons while preserving exact counts for retained repeated keys. It can therefore change low-coverage reconstruction relative to an unfiltered graph. Yield and error counts accompany every timing. The alternative no-Bloom c=1 point helps distinguish singleton pruning from filter implementation cost in the c=0 datasets.

`format_graph_bloom.c` changes the format caller to invoke the existing `count_files_filtered` and `c4x_set_bloom` APIs. The saved counter, canonical graph, bounded traversal and Bloom implementation remain the same. The original no-Bloom caller and its complete measurements are retained. Source and binary hashes document both callers. Separate smoke tests verify no-Bloom, B8, B16, F32, F34 and B12/F34 on 100 known noise-free designs of each format.

## Conventional comparator and quality

This is a reference-free implementation comparator, not a reproduction of the original full LSH or BMA decoders. Exact observed prefix keys are 16 nt for Fountain/HT176 and 14 nt for Organick. Each Organick molecule is normalized to one canonical orientation, avoiding duplicate forward/reverse alignments. HT176 uses the same 15-nt primer suffix and fixed numeric index range as graph seed construction; unknown designed-address membership is never looked up.

Eligible clusters have at least three observations. One MUSCLE alignment uses at most 15 reads. Fountain/Organick use deterministic uniform reservoir sampling. For fragmented HT176 reads, selection prefers the longest observed coverage of the known 164-nt region and uses random tie-breaking. Terminal alignment gaps are missing observations; internal gaps retain the 0.4 vote weight. No ground-truth sequence determines selection, voting or conflict resolution. A dedicated test reconstructs all 100 designs from one complete observation plus 19 short prefix fragments per design.

Candidates use the same format validation as graph outputs: Fountain RS(38,36) in detect-only mode, the globally calibrated Organick rotating-ternary predicate already described in the manuscript, and the HT176 legacy two-base mapping with CRC16. Conflicting valid candidates for a prefix are rejected without consulting designs. Exact-reference scoring follows the timer. Organick does not have an implemented outer decoder here. Fountain code validity and complete-file checksum recovery are separate endpoints.

## Timing and resource interpretation

The timer spans subprocess startup, input parsing, clustering or graph counting, seeding, MSA/voting or path search, code validation, reference-free conflict rejection and sorted output writing. Query-filter construction and both counting passes are included. Downloads, shared preparation, tuning, truth scoring and outer decoding are excluded. Shared merge-and-conversion preparation was measured once: Fountain 46.28 s and HT176 142.29 s. The Fountain graph view additionally prepends a synthetic constant boundary anchor and appends one synthetic terminal lookahead; no payload bases are supplied or changed. Its view-creation cost is excluded from reconstruction and measured in a separate post-matrix replay, recorded in annotation_preparation.json. Other methods consume the same retained payload molecules in the plain view. Add these costs when budgeting from original paired FASTQ; the plotted speed ratios concern reconstruction after common preprocessing.

All measured processes run serially on the same Apple M3 Max workstation. Scheduled method order rotates by dataset/repeat. Completed unchanged comparator measurements from the documented earlier block are retained; no partial or interrupted result enters the matrix. Graph counting uses eight insertion workers plus reading/extraction pipeline stages; search is single-threaded. The baseline uses eight single-threaded alignment workers plus coordinator overhead. Total active threads/core use are not constrained to equality. Numerical-library threading is capped at one. Background applications remain present; no CPU affinity, dedicated idle host or forced cold-cache regime is claimed. Pre-run process snapshots and complete time(1) stderr are retained.

Speed ratios divide mean conventional-baseline reconstruction time by mean configuration time. They are operating-point comparisons, not universal, equal-recovery or complete-outer-decoder speedups. Standard deviations describe computational repeat variability. Repeated output hashes verify deterministic quality counts.

Peak RSS in the memory figure is the maximum native graph-executable RSS across three repeats, in decimal GB. It includes bitmap/query-filter construction phases. The time(1) maximum for an MSA process is not summed process-tree memory, so it is excluded from cross-configuration memory plots.

## Reproduction and evidence

- `prepare.py`: shared NGmerge/FASTA preparation; original paths can be edited for another host.
- `provenance.py`: raw/processed/reference sizes and hashes.
- `run_benchmark.py --pilot --only song_ht176`: cutoff pilot (unchanged public-file settings retained).
- `smoke_test.py`, `fragment_smoke_test.py`, `bloom_smoke_test.py`: adapter/terminal-gap/Bloom contract tests, not real-data performance estimates.
- `run_bloom_matrix.py`: complete serial timing matrix, resuming only from complete result JSON.
- `validate_matrix.py`: completion, hashes, query-filter invariance, counting-mode consistency, and untimed Fountain file/payload checks.
- `summarize_matrix.py`: source CSVs and vector/raster Figures 10–11.
- `update_matrix_manuscript.py`: OOXML-preserving manuscript update after measurements and validation.

Large reads/FASTA files are excluded from the portable archive. Canonical outputs can be deduplicated by SHA256 while all 81 repeat records, exact commands, stage logs and output hashes are retained. Development/rejected runs are clearly separated and excluded. Production source and manuscript history were not committed.
