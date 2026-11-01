# File 1: complete decoding experiment report

Version: 2026-10-18. Includes all recorded stages in this task. Numerical tables are generated from local JSON/TSV records; both languages use the same data.

## 1. Outcome and scope

File 1 was decoded from actual strand decisions generated from sequencing reads. The best strict graph output contains 16,125 correct strands and 0 wrong strands. The recovered ZIP passes byte equality and internal CRC validation. No design strands were injected.

DBGPS2 here means an experimental adapter using this repository's exact canonical k-mer counter, graph queries and dynamic coverage-ratio rule, with an internal-index bidirectional traversal and GF(64) validation. It is not direct support by the stock CLI. Reference payloads do not seed graph traversal. Reference assistance in the later File 3 read-routing step is disclosed separately.

## 2. Inputs, integrity and preprocessing

| Item | Value |
| --- | --- |
| FASTQ | I16_S2_R1_001.fastq.gz |
| Figshare file ID | 22305708 |
| Compressed bytes | 994447128 |
| MD5 | 6d48dee09e7cbe3f11bc749c9a7839c6 |
| SHA-256 | 25b70f73fcb4582666541f5379a74dc86d38e5e4ec690d8fad002986b03dad48 |
| Raw reads | 29303855 |
| Retained reads | 14299688 |
| Design strands | 16383 |
| Outer blocks | 1 |

The complete compressed input passed its published MD5. Baseline graph and matched MSA runs use every read of original length 55-70 nt; reads longer than 60 nt are capped at the first 60 nt, and shorter retained reads keep their length. These are not subsampled experiments. The sample_reads metadata describes an additional pilot file that is not used for the reported full-data results.

Physical File 1 designs and the original ZIP come from the authors' repository. The later public design copy agrees with the primary design. Re-encoding reproduces all 16,383 physical strands exactly.

## 3. Encoding and outer-code criterion

The physical layout is 42-nt payload + 12-nt scrambled index + 6-nt inner parity, totaling 60 nt. A/C/G/T map to 0/1/2/3. The 24-bit index is XOR-scrambled with 1071068. Inner RS(20,18) over GF(2^6), polynomial 91, corrects one 6-bit symbol. It is not a general indel decoder: traversal produces full-length candidates before symbol correction.

Each block contains six independent RS(16383,10977) codewords over GF(2^14), with 5,406 parity symbols. A strand carries one 14-bit symbol in each lane. Every lane must satisfy 2E_j + S <= 5406, where S is missing strands and E_j is wrong symbols in that lane. With no errors, 10,977 distinct indices suffice. Wrong-strand counts give a conservative bound; later results use actual lane errors. Blocks are independent, with no additional cross-block parity.

## 4. Traversal, metric definitions and controls

Seeds are observed graph k-mers containing a valid 12-nt index. Traversal anchors at internal-index end position 54, extends left to the strand start, then right through the final 6 nt. At most 1,050 paths remain active per generation; overflow erases the entire source-index search. Strict output requires a unique full-length inner-code-valid strand. Canonical counts combine reverse complements and saturate at 16,383.

Cutoff c means count > c, not >= c. For previous coverage p, next coverage q must be within [max(c+1, floor(p/r)), min(16383, floor(p*r))]. All measured grids use c>=1; the two-pass singleton filter preserves eligible counts >=2. Changing c alone to zero would not retain singletons because seed enumeration also requires count>1.

Four outcomes are kept separate: (1) exact strands directly observed in input reads; (2) indices with at least one exact full 60-nt candidate; (3) correct, wrong and erased strands after an actual one-sequence-per-index decision; and (4) complete-file bytes and ZIP validation. Cross-parameter candidate unions require truth labels and are not an implemented reference-free selector.

Stock-CLI diagnostic runs emit zero strands on both real and perfect-reference inputs because their prefix-index, primer and CRC/byte-RS semantics do not match this format. Compatibility diagnosis established the historical libc++ uniform sampler and spectral-message outer decoding. A perfect-reference k=21 positive control recovers all 16,383 strands and the original ZIP. The best real graph output is byte-identical when no evaluation FASTA is supplied.

## 5. Pooled-read MSA and graph attempts

Matched-input MSA uses a 16-base prefix, minimum cluster size 5, MUSCLE 3.8.31 and gap weight 0.4. Large clusters produce five subsets of at most 15 reads with fixed seed 20261017. This is a prefix-clustering baseline, not a full reproduction of the paper's LSH procedure. Recorded timings reflect concurrent workloads and differing workers, not controlled speed comparisons.

| Method | Correct indices | Wrong outputs | Interpretation |
| --- | --- | --- | --- |
| Direct exact input | 4439 | - | 27.0952% |
| MSA candidate coverage | 10779 | - | 65.7938% |
| MSA strict decisions | 10778 | 4 | GF(64) check |
| MSA corrected decisions | 13179 | 779 | One-symbol correction |
| Graph strict decisions | 16125 | 0 | k19/c1/r2.55 |

MSA processes 215,090 clusters and generates 350,470 candidates, of which 346,253 are 60 nt. Recorded total processing time is 281.007 s, not a controlled cross-algorithm benchmark. Inner correction can improve correct decisions while also introducing miscorrections, so wrong outputs are reported alongside gains.

An additional run follows the released notebook: prefix 14, minimum cluster 4, gap weight 0.5, alignment of complete 55-70-nt reads before capping the voted candidate at 60 nt. Its preprocessing differs from the matched baseline and is reported separately.

| Notebook variant | Correct | Wrong |
| --- | --- | --- |
| Candidate coverage | 11938 | - |
| Strict decisions | 11934 | 3 |
| Corrected decisions | 13602 | 706 |

The notebook strict output also independently decodes File 1 to byte-identical ZIP content with valid CRCs; recorded outer-decode time is 2.243 s.

The final pooled strict grid contains 79 settings. The best strict output at each k is listed below. All individual settings and saved commands, exit codes and timings are included in evidence/File1/graph_parameters.csv and run_ledger.json.

| k | Settings | Best c | Best r | Correct | Wrong | Overflow indices |
| --- | --- | --- | --- | --- | --- | --- |
| 12 | 4 | 1 | 10.0 | 0 | 0 | 16383 |
| 13 | 9 | 199 | 2.55 | 65 | 1 | 7480 |
| 15 | 12 | 19 | 2.55 | 11091 | 1 | 3440 |
| 17 | 12 | 3 | 2.55 | 16047 | 0 | 24 |
| 19 | 10 | 1 | 2.55 | 16125 | 0 | 0 |
| 21 | 12 | 1 | 2.55 | 16006 | 0 | 0 |
| 25 | 10 | 1 | 2.55 | 15495 | 0 | 0 |
| 31 | 10 | 1 | 2.55 | 13983 | 0 | 0 |

Smaller k increases local support but also cross-strand sharing and noisy branches; larger k loses low-support paths. True-path support alone therefore does not establish actual recoverability. Best parameters are selected retrospectively against this dataset's truth, without a held-out validation set.

## 6. Read-partition status and per-block verification

File 1 has one outer-code block, so there is no cross-block read-partition comparison. Reference-alignment read routing was not run for this file. All 79 graph settings use the retained pooled input; File 3 partitioning gains are not transferred to File 1.

| Block (1-based) | Correct | Wrong | Erased | Conservative 2E+S | Pass |
| --- | --- | --- | --- | --- | --- |
| 1 | 16125 | 0 | 258 | 258 | True |

All 1 blocks pass. In the scan, 30 settings satisfy the conservative bound in every block. This does not mean that each was fully decoded; actual file-decoding attempts are identified by the decode records.

## 7. Complete-file and implementation validation

| Validation | Result |
| --- | --- |
| Outer decoder exit | 0 |
| Decoded bytes | 99103 |
| SHA-256 | 05295bc085b243e16aafa18a63a7b5fad944037f68667a7ac60bff691ac15945 |
| Byte equality | True |
| ZIP CRC error | None |
| Outer seconds | 0.291609 |

rs_driver_fast uses the authors' spectral outer-code algorithm and historical randomization with an equivalent GF(2^14) log-table acceleration. Evidence includes 232,072 finite-field checks and original-versus-accelerated File 1 byte comparison. File success is established from the decoded ZIP, not fragment recall alone.

Full output: artifacts/File1_sequencing_recovered.zip. File 1 is compared to the authors' original ZIP; File 2 is compared to the CRC-valid ZIP decoded from the complete physical reference.

## Appendix A. Exact grid and execution records

The table lists measured c values for each k/r present in final pooled metrics. Unlisted combinations are not implied to have been tested. Evidence contains individual correctness, errors, overflow, ambiguity, block budgets, commands and timings. Earlier repeats may be superseded by later runs; ledgers retain history, while final metrics count each setting once.

| k | r | Measured c values |
| --- | --- | --- |
| 12 | 2.55 | 1, 3 |
| 12 | 10.0 | 1, 3 |
| 13 | 2.55 | 1, 3, 9, 19, 49, 99, 199 |
| 13 | 10.0 | 1, 3 |
| 15 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 15 | 10.0 | 1, 3 |
| 17 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 17 | 10.0 | 1, 3 |
| 19 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 21 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 21 | 10.0 | 1, 3 |
| 25 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |
| 31 | 2.55 | 1, 3, 4, 5, 6, 9, 19, 49, 99, 199 |

## Appendix B. Reproduction entry points and script roles

Standalone code archive: File1_code.zip. Preserve its layout and run from its root. It contains the recorded DBGPS2 core, adapters, portable RS-header snapshot, MUSCLE source archive, pinned-version metadata and small validation fixtures. Raw reads are excluded; the explicit fetch stage downloads the matching complete FASTQ and verifies MD5.

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

Primary local evidence is under evidence/File1/: original summary.json, graph_metrics.json, graph_parameters.csv, run_ledger.json, verified-input records and decode validation. MANIFEST.sha256 records delivered-file digests; EVIDENCE_MAP.json links each evidence item to its original workspace path.

[Article: DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3). [Primary repository](https://github.com/MLI-lab/noisy_dna_data_storage), pinned to 41422be448b3c1da9041c3a2967d2754d2ad8b1a. [RS repository](https://github.com/reinhardh/dna_rs_coding), checkout 455e1a5182bf64281f54d434a4ed5ca456c2373e, historical headers from 8e7e3b5. [Later designs](https://github.com/fml-ethz/dt4dds-challenges_notebooks), pinned to 861d3b95b52d1b121e9d85040f20c8ecf0cda2d4.
