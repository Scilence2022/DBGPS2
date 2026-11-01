# File 2: complete decoding experiment report

Version: 2026-10-18. Includes all recorded stages in this task. Numerical tables are generated from local JSON/TSV records; both languages use the same data.

## 1. Outcome and scope

File 2 was decoded from actual strand decisions generated from sequencing reads. The best strict graph output contains 43,472 correct strands and 3 wrong strands. The recovered ZIP passes byte equality and internal CRC validation. No design strands were injected.

DBGPS2 here means an experimental adapter using this repository's exact canonical k-mer counter, graph queries and dynamic coverage-ratio rule, with an internal-index bidirectional traversal and GF(64) validation. It is not direct support by the stock CLI. Reference payloads do not seed graph traversal. Reference assistance in the later File 3 read-routing step is disclosed separately.

## 2. Inputs, integrity and preprocessing

| Item | Value |
| --- | --- |
| FASTQ | I18_S3_R1_001.fastq.gz |
| Figshare file ID | 22305990 |
| Compressed bytes | 3150366942 |
| MD5 | 3e76a04053fde5c0ee43e81d2c50c807 |
| SHA-256 | f7aee402208bd6d6fcc31d6d303d178ba787f8b12f4be740f28d87616ac8ac31 |
| Raw reads | 100706815 |
| Retained reads | 33404713 |
| Design strands | 49149 |
| Outer blocks | 3 |

The complete compressed input passed its published MD5. Baseline graph and matched MSA runs use every read of original length 55-70 nt; reads longer than 60 nt are capped at the first 60 nt, and shorter retained reads keep their length. These are not subsampled experiments. The sample_reads metadata describes an additional pilot file that is not used for the reported full-data results.

Later File 2/3 design deposits do not use the physical strand orientation and index layout. prepare_references.py reverses base order, not reverse complements, preserves the 42-nt payload, restores the scrambled internal index and regenerates inner parity. This transformation is validated on every primary File 1 strand, and the transformed reference for this file yields a ZIP with valid CRCs.

A reference-format check on the first 1,000,000 raw reads finds 88 full exact read matches. That subset is only a layout validation; reconstruction uses the complete retained input. The byte-comparison baseline is decoded from the complete validated design, not assumed to have been recovered from incomplete sequencing output.

## 3. Encoding and outer-code criterion

The physical layout is 42-nt payload + 12-nt scrambled index + 6-nt inner parity, totaling 60 nt. A/C/G/T map to 0/1/2/3. The 24-bit index is XOR-scrambled with 1071068. Inner RS(20,18) over GF(2^6), polynomial 91, corrects one 6-bit symbol. It is not a general indel decoder: traversal produces full-length candidates before symbol correction.

Each block contains six independent RS(16383,10977) codewords over GF(2^14), with 5,406 parity symbols. A strand carries one 14-bit symbol in each lane. Every lane must satisfy 2E_j + S <= 5406, where S is missing strands and E_j is wrong symbols in that lane. With no errors, 10,977 distinct indices suffice. Wrong-strand counts give a conservative bound; later results use actual lane errors. Blocks are independent, with no additional cross-block parity.

## 4. Traversal, metric definitions and controls

Seeds are observed graph k-mers containing a valid 12-nt index. Traversal anchors at internal-index end position 54, extends left to the strand start, then right through the final 6 nt. At most 1,050 paths remain active per generation; overflow erases the entire source-index search. Strict output requires a unique full-length inner-code-valid strand. Canonical counts combine reverse complements and saturate at 16,383.

Cutoff c means count > c, not >= c. For previous coverage p, next coverage q must be within [max(c+1, floor(p/r)), min(16383, floor(p*r))]. All measured grids use c>=1; the two-pass singleton filter preserves eligible counts >=2. Changing c alone to zero would not retain singletons because seed enumeration also requires count>1.

Four outcomes are kept separate: (1) exact strands directly observed in input reads; (2) indices with at least one exact full 60-nt candidate; (3) correct, wrong and erased strands after an actual one-sequence-per-index decision; and (4) complete-file bytes and ZIP validation. Cross-parameter candidate unions require truth labels and are not an implemented reference-free selector.

Early File 2 diagnostics used the later deposited base order directly and also explored a prefix-index interpretation. Those artifacts remain under orientation_diagnostic and are excluded from the final 35-setting grid; they are not final recovery rates. Physical-layout checks, exact raw-read matches and the reference ZIP jointly validate the corrected conversion.

## 5. Pooled-read MSA and graph attempts

Matched-input MSA uses a 16-base prefix, minimum cluster size 5, MUSCLE 3.8.31 and gap weight 0.4. Large clusters produce five subsets of at most 15 reads with fixed seed 20261017. This is a prefix-clustering baseline, not a full reproduction of the paper's LSH procedure. Recorded timings reflect concurrent workloads and differing workers, not controlled speed comparisons.

| Method | Correct indices | Wrong outputs | Interpretation |
| --- | --- | --- | --- |
| Direct exact input | 2802 | - | 5.7010% |
| MSA candidate coverage | 11270 | - | 22.9303% |
| MSA strict decisions | 11264 | 13 | GF(64) check |
| MSA corrected decisions | 18287 | 8001 | One-symbol correction |
| Graph strict decisions | 43472 | 3 | k21/c1/r2.55 |

MSA processes 388,269 clusters and generates 617,541 candidates, of which 601,364 are 60 nt. Recorded total processing time is 1085.370 s, not a controlled cross-algorithm benchmark. Inner correction can improve correct decisions while also introducing miscorrections, so wrong outputs are reported alongside gains.

The final pooled strict grid contains 35 settings. The best strict output at each k is listed below. All individual settings and saved commands, exit codes and timings are included in evidence/File2/graph_parameters.csv and run_ledger.json.

| k | Settings | Best c | Best r | Correct | Wrong | Overflow indices |
| --- | --- | --- | --- | --- | --- | --- |
| 13 | 5 | 1 | 2.55 | 0 | 0 | 49117 |
| 15 | 5 | 19 | 2.55 | 14018 | 20 | 23828 |
| 17 | 5 | 9 | 2.55 | 36620 | 2 | 444 |
| 19 | 5 | 1 | 2.55 | 43007 | 7 | 2351 |
| 21 | 5 | 1 | 2.55 | 43472 | 3 | 0 |
| 25 | 5 | 1 | 2.55 | 37673 | 3 | 0 |
| 31 | 5 | 1 | 2.55 | 24086 | 4 | 0 |

Smaller k increases local support but also cross-strand sharing and noisy branches; larger k loses low-support paths. True-path support alone therefore does not establish actual recoverability. Best parameters are selected retrospectively against this dataset's truth, without a held-out validation set.

## 6. Read-partition status and per-block verification

File 2 has three outer-code blocks. Design-alignment read partitioning was not performed for File 2, so no post-partition measurement is available. The pooled result already satisfies the bound in all three blocks and recovers the full ZIP. The table groups decoded indices by block; it does not imply that reads were partitioned before graph decoding.

| Block (1-based) | Correct | Wrong | Erased | Conservative 2E+S | Pass |
| --- | --- | --- | --- | --- | --- |
| 1 | 14429 | 2 | 1952 | 1956 | True |
| 2 | 15443 | 1 | 939 | 941 | True |
| 3 | 13600 | 0 | 2783 | 2783 | True |

All 3 blocks pass. In the scan, 6 settings satisfy the conservative bound in every block. This does not mean that each was fully decoded; actual file-decoding attempts are identified by the decode records.

## 7. Complete-file and implementation validation

| Validation | Result |
| --- | --- |
| Outer decoder exit | 0 |
| Decoded bytes | 320949 |
| SHA-256 | 4664730829f1a61003a18dd276db736bc13cc2de647bff012cdb1493895b0cd5 |
| Byte equality | True |
| ZIP CRC error | None |
| Outer seconds | 2.827647 |

rs_driver_fast uses the authors' spectral outer-code algorithm and historical randomization with an equivalent GF(2^14) log-table acceleration. Evidence includes 232,072 finite-field checks and original-versus-accelerated File 1 byte comparison. File success is established from the decoded ZIP, not fragment recall alone.

Full output: artifacts/File2_sequencing_recovered.zip. File 1 is compared to the authors' original ZIP; File 2 is compared to the CRC-valid ZIP decoded from the complete physical reference.

## Appendix A. Exact grid and execution records

The table lists measured c values for each k/r present in final pooled metrics. Unlisted combinations are not implied to have been tested. Evidence contains individual correctness, errors, overflow, ambiguity, block budgets, commands and timings. Earlier repeats may be superseded by later runs; ledgers retain history, while final metrics count each setting once.

| k | r | Measured c values |
| --- | --- | --- |
| 13 | 2.55 | 1, 3, 9, 19, 49 |
| 15 | 2.55 | 1, 3, 9, 19, 49 |
| 17 | 2.55 | 1, 3, 9, 19, 49 |
| 19 | 2.55 | 1, 3, 9, 19, 49 |
| 21 | 2.55 | 1, 3, 9, 19, 49 |
| 25 | 2.55 | 1, 3, 9, 19, 49 |
| 31 | 2.55 | 1, 3, 9, 19, 49 |

## Appendix B. Reproduction entry points and script roles

Standalone code archive: File2_code.zip. Preserve its layout and run from its root. It contains the recorded DBGPS2 core, adapters, portable RS-header snapshot, MUSCLE source archive, pinned-version metadata and small validation fixtures. Raw reads are excluded; the explicit fetch stage downloads the matching complete FASTQ and verifies MD5.

```sh
python3 tools/reproduce.py build
python3 tools/reproduce.py smoke
python3 tools/reproduce.py fetch
python3 tools/reproduce.py prepare
python3 tools/reproduce.py graph
python3 tools/reproduce.py msa
```

```sh
python3 tools/reproduce.py outer
```

| Script | Role |
| --- | --- |
| download.py / fetch_sources.py | Verified downloads and source pins |
| prepare_reads.c / prepare_references.py | Length preprocessing and reference transcoding |
| graph_adapter.c / run_graphs.py | Strict traversal and parameter grid |
| msa_baseline.py / inner_select.cpp | MSA and common inner-code decisions |
| rs_driver.cpp / fast_gf14.hpp / legacy_random.hpp | Outer decode, field acceleration and legacy sampler |

## Appendix C. Limitations, sources and evidence map

All best settings are exploratory choices on this dataset, with no held-out validation or multi-seed confidence intervals. MSA is not a complete LSH reproduction. Coverage, candidate unions, actual decisions and file recovery are distinct outcomes. Filtering, length handling, partitioning and added truth information are disclosed separately. Unattempted, interrupted and diagnostic failures are not represented as completed measurements.

Primary local evidence is under evidence/File2/: original summary.json, graph_metrics.json, graph_parameters.csv, run_ledger.json, verified-input records and decode validation. MANIFEST.sha256 records delivered-file digests; EVIDENCE_MAP.json links each evidence item to its original workspace path.

[Article: DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3). [Primary repository](https://github.com/MLI-lab/noisy_dna_data_storage), pinned to 41422be448b3c1da9041c3a2967d2754d2ad8b1a. [RS repository](https://github.com/reinhardh/dna_rs_coding), checkout 455e1a5182bf64281f54d434a4ed5ca456c2373e, historical headers from 8e7e3b5. [Later designs](https://github.com/fml-ethz/dt4dds-challenges_notebooks), pinned to 861d3b95b52d1b121e9d85040f20c8ecf0cda2d4.
