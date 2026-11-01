# File 3: complete decoding experiment report

Version: 2026-10-18. Includes all recorded stages in this task. Numerical tables are generated from local JSON/TSV records; both languages use the same data.

## 1. Outcome and scope

No complete File 3 has been recovered from the measured sequencing outputs alone. After reference-assisted read partitioning, actual outputs for blocks 1-11 satisfy the outer-code bound; block 12 does not. Two explicitly design-supplemented simulations do recover the complete file, and are reported separately.

DBGPS2 here means an experimental adapter using this repository's exact canonical k-mer counter, graph queries and dynamic coverage-ratio rule, with an internal-index bidirectional traversal and GF(64) validation. It is not direct support by the stock CLI. Reference payloads do not seed graph traversal. Reference assistance in the later File 3 read-routing step is disclosed separately.

## 2. Inputs, integrity and preprocessing

| Item | Value |
| --- | --- |
| FASTQ | File3.fq.gz |
| Figshare file ID | 49503285 |
| Compressed bytes | 5723104027 |
| MD5 | b08fb297825ee97adee126c900bd50c6 |
| SHA-256 | 98ca69f7532f8ccbd7d57644dbe4830dec52e0b173d74d5485f37f521b697303 |
| Raw reads | 195619515 |
| Retained reads | 92337721 |
| Design strands | 196596 |
| Outer blocks | 12 |

The complete compressed input passed its published MD5. Baseline graph and matched MSA runs use every read of original length 55-70 nt; reads longer than 60 nt are capped at the first 60 nt, and shorter retained reads keep their length. These are not subsampled experiments. The sample_reads metadata describes an additional pilot file that is not used for the reported full-data results.

Later File 2/3 design deposits do not use the physical strand orientation and index layout. prepare_references.py reverses base order, not reverse complements, preserves the 42-nt payload, restores the scrambled internal index and regenerates inner parity. This transformation is validated on every primary File 1 strand, and the transformed reference for this file yields a ZIP with valid CRCs.

A reference-format check on the first 1,000,000 raw reads finds 108 full exact read matches. That subset is only a layout validation; reconstruction uses the complete retained input. The byte-comparison baseline is decoded from the complete validated design, not assumed to have been recovered from incomplete sequencing output.

## 3. Encoding and outer-code criterion

The physical layout is 42-nt payload + 12-nt scrambled index + 6-nt inner parity, totaling 60 nt. A/C/G/T map to 0/1/2/3. The 24-bit index is XOR-scrambled with 1071068. Inner RS(20,18) over GF(2^6), polynomial 91, corrects one 6-bit symbol. It is not a general indel decoder: traversal produces full-length candidates before symbol correction.

Each block contains six independent RS(16383,10977) codewords over GF(2^14), with 5,406 parity symbols. A strand carries one 14-bit symbol in each lane. Every lane must satisfy 2E_j + S <= 5406, where S is missing strands and E_j is wrong symbols in that lane. With no errors, 10,977 distinct indices suffice. Wrong-strand counts give a conservative bound; later results use actual lane errors. Blocks are independent, with no additional cross-block parity.

## 4. Traversal, metric definitions and controls

Seeds are observed graph k-mers containing a valid 12-nt index. Traversal anchors at internal-index end position 54, extends left to the strand start, then right through the final 6 nt. At most 1,050 paths remain active per generation; overflow erases the entire source-index search. Strict output requires a unique full-length inner-code-valid strand. Canonical counts combine reverse complements and saturate at 16,383.

Cutoff c means count > c, not >= c. For previous coverage p, next coverage q must be within [max(c+1, floor(p/r)), min(16383, floor(p*r))]. All measured grids use c>=1; the two-pass singleton filter preserves eligible counts >=2. Changing c alone to zero would not retain singletons because seed enumeration also requires count>1.

Four outcomes are kept separate: (1) exact strands directly observed in input reads; (2) indices with at least one exact full 60-nt candidate; (3) correct, wrong and erased strands after an actual one-sequence-per-index decision; and (4) complete-file bytes and ZIP validation. Cross-parameter candidate unions require truth labels and are not an implemented reference-free selector.

## 5. Pooled-read MSA and graph attempts

Matched-input MSA uses a 16-base prefix, minimum cluster size 5, MUSCLE 3.8.31 and gap weight 0.4. Large clusters produce five subsets of at most 15 reads with fixed seed 20261017. This is a prefix-clustering baseline, not a full reproduction of the paper's LSH procedure. Recorded timings reflect concurrent workloads and differing workers, not controlled speed comparisons.

| Method | Correct indices | Wrong outputs | Interpretation |
| --- | --- | --- | --- |
| Direct exact input | 8199 | - | 4.1705% |
| MSA candidate coverage | 41564 | - | 21.1418% |
| MSA strict decisions | 41560 | 19 | GF(64) check |
| MSA corrected decisions | 52893 | 14464 | One-symbol correction |
| Graph strict decisions | 137662 | 143 | k19/c1/r3.5 |

MSA processes 365,025 clusters and generates 654,277 candidates, of which 641,449 are 60 nt. Recorded total processing time is 895.230 s, not a controlled cross-algorithm benchmark. Inner correction can improve correct decisions while also introducing miscorrections, so wrong outputs are reported alongside gains.

The final pooled strict grid contains 179 settings. The best strict output at each k is listed below. All individual settings and saved commands, exit codes and timings are included in evidence/File3/graph_parameters.csv and run_ledger.json.

| k | Settings | Best c | Best r | Correct | Wrong | Overflow indices |
| --- | --- | --- | --- | --- | --- | --- |
| 13 | 14 | 1 | 2.55 | 0 | 0 | 196584 |
| 15 | 25 | 19 | 1.8 | 9061 | 24 | 128951 |
| 17 | 25 | 9 | 3.5 | 69758 | 13 | 5444 |
| 18 | 25 | 3 | 3.0 | 115984 | 18 | 8839 |
| 19 | 34 | 1 | 3.5 | 137662 | 143 | 9150 |
| 20 | 21 | 1 | 3.5 | 136445 | 61 | 67 |
| 21 | 25 | 1 | 3.5 | 126757 | 58 | 69 |
| 25 | 5 | 1 | 2.55 | 83651 | 39 | 1 |
| 31 | 5 | 1 | 2.55 | 35522 | 37 | 21 |

Smaller k increases local support but also cross-strand sharing and noisy branches; larger k loses low-support paths. True-path support alone therefore does not establish actual recoverability. Best parameters are selected retrospectively against this dataset's truth, without a held-out validation set.

## 6. Before partitioning: smaller k and ratio expansion

The original 35 settings at r=2.55 peak at k19/c1: 135,575 correct and 111 wrong strands. A further 144 settings expand the grid to 179, including k18/k20, smaller k, and ratios from 1.2 to 10. The best strict correct count is k19/c1/r3.5: 137,662 correct, 143 wrong, with 8/12 blocks in bound. Maximizing passing blocks instead selects k20/c1/r3.5: 136,445 correct, 61 wrong, and 9/12 blocks in bound.

| Block | Correct | Wrong | Erased | Conservative budget | Pass |
| --- | --- | --- | --- | --- | --- |
| 1 | 10694 | 17 | 5672 | 5706 | False |
| 2 | 12972 | 6 | 3405 | 3417 | True |
| 3 | 13265 | 7 | 3111 | 3125 | True |
| 4 | 12636 | 10 | 3737 | 3757 | True |
| 5 | 13368 | 6 | 3009 | 3021 | True |
| 6 | 10425 | 18 | 5940 | 5976 | False |
| 7 | 12895 | 9 | 3479 | 3497 | True |
| 8 | 12023 | 12 | 4348 | 4372 | True |
| 9 | 11169 | 11 | 5203 | 5225 | True |
| 10 | 11198 | 12 | 5173 | 5197 | True |
| 11 | 10824 | 13 | 5546 | 5572 | False |
| 12 | 6193 | 22 | 10168 | 10212 | False |

Block 12 has only 6,193 correct and 22 wrong strands in that best strict result, with 10,168 erasures and a conservative budget of 10,212. Although k13/k15 increase true-k-mer support, extensive branching causes overflow. Graph support is not equivalent to emitted reconstruction.

## 7. Candidate-level inner RS and the complete parameter union

One-symbol RS correction is applied to every complete, non-overflow graph path, not reapplied to strict outputs that already pass parity. Repaired indices can be reassigned; out-of-domain indices are rejected. Eleven actual rules are evaluated: strict, corrected_unique, corrected_vote1/2/3/5/10/20, and pooled_weight1/4/16. Votes count graph paths, not reads. Priority rules prefer a unique zero-correction word; other decisions use uniqueness or vote thresholds, with ties erased. Pooled rules choose a unique maximum using combined weighted zero/one-correction votes.

Initially a complete-path-count upper bound excludes 119 settings from improving the attained single-setting maximum, leaving 60 reruns and 660 decision outputs. This pruning is valid for a single-setting maximum, not for a union. All 119 omitted settings are therefore later run, completing the original 179. Strict outputs are checked byte-for-byte against the earlier version.

| Metric | Result |
| --- | --- |
| Best single candidate coverage | 140797 |
| Best actual unique correct | 138674 |
| Wrong at that actual output | 17144 |
| 179-setting correct candidate union | 151202 |
| Zero-correction union | 148835 |
| RS required across all settings | 2367 |
| All-policy output correct union | 150526 |
| Block12 candidate union | 6927 |

At k19/c1/r3.5, 3,361,623 complete paths comprise 139,101 zero-correction valid paths, 1,000,579 one-symbol-corrected paths, 2,161,363 uncorrectable paths and 60,580 out-of-domain paths. Correct candidate coverage is 140,797, including 2,496 indices requiring correction at this setting. The best actual correct count is k19/c1/r3/pooled_weight16: 138,674 correct but 17,144 wrong, passing only 5/12 blocks. Correct count alone is not a file-recovery criterion.

| Block | 179-setting candidate union | Missing | Oracle erasure bound |
| --- | --- | --- | --- |
| 1 | 11383 | 5000 | True |
| 2 | 13619 | 2764 | True |
| 3 | 14032 | 2351 | True |
| 4 | 13340 | 3043 | True |
| 5 | 14156 | 2227 | True |
| 6 | 14091 | 2292 | True |
| 7 | 13951 | 2432 | True |
| 8 | 12805 | 3578 | True |
| 9 | 12896 | 3487 | True |
| 10 | 12371 | 4012 | True |
| 11 | 11631 | 4752 | True |
| 12 | 6927 | 9456 | False |

The union counts each repaired decoded index once if its full exact 60-nt truth occurs under any setting: 151,202/196,596 (76.9100%). Truth-selected unions for blocks 1-11 satisfy an erasure-only criterion; block 12 still misses 9,456 strands. This is candidate coverage, not a reference-free cross-parameter selection algorithm or complete File 3 recovery.

## 8. DBGPS-analyzer diagnosis of block 12

The user-specified analyzer is compiled from its recorded source without changes. Thirteen full count runs cover nine k values on processed reads, k17/k19/k20 on original full-length reads, and original reads capped at 60 nt at k19. DBGPS-links analyzes the entire design and each of the twelve blocks at every k. Targets remain full 60-nt strands with no primer trimming.

| k19 input | Block12 all k-mers >1 | Median strand mean coverage |
| --- | --- | --- |
| Processed 55-70/cap60 | 6277 | 4.12 |
| All raw/cap60 | 7768 | 5.64 |
| All raw/full length | 7901 | 5.71 |

Among 9,925 block-12 failures after RS at k19/c1/r3.5, 7,430 have an absent true k-mer, 2,429 have only singleton minimum support, 28 fail the ratio rule, and 38 have an eligible true path but overflow. The full-library k19 distinct-noise/target-key ratio Kn=27.516543 is not a base error rate. Block12 has 50.85% GC and three within-block design cross-links at k19 (other blocks: 1-7), without an obvious coarse design anomaly. Pooled canonical coverage is not target-assigned read depth.

Per-index minimum coverage agrees with the original counter at all nine k values. Directed integer-rounded feasibility plus recorded overflow states exactly predicts zero-correction correct-candidate flags across all 179 settings. This validates the diagnosis, not a claim that all feasible truth paths were emitted.

## 9. Design-alignment read routing and the same-read re-pooling control

minimap2 2.31-r1302 aligns all 195,619,515 raw reads to the complete design with -a -x sr -k11 -w2 -n2 -m15 -s20 -A2 -B3 -O4,12 -E2,1 -N5 -p0.8 --secondary=yes --eqx -t12. Full payloads assist assignment, so this is reference-assisted. Routing changes no read bases or gaps. Reverse-strand SAM primaries are restored to original sequence orientation and quality order.

Broad assignment requires >=20 matches, >=30 aligned query bases, >=20 payload-overlap bases, AS>=25 and >=60% match fraction. Every reported qualifying alignment scoring at least 80% of the best qualifying hit must name the same block. Strict additionally requires primary MAPQ>=10, >=24 matches, >=35 query bases, >=25 payload bases, AS>=40 and >=70% match fraction. At most five secondary hits are reported; unreported alternatives may still exist.

Broad pools contain 44,581,669 reads, with 15,770 cross-block ambiguous and 151,022,076 unassigned reads. Twelve broad pools plus ambiguous and unassigned reads form a complete disjoint partition; strict/inclusive variants overlap them. Full read-name/sequence/quality SHA-256 multiset totals, sums and XORs agree with the source, validating complete record conservation rather than only a pilot subset.

| Block | Broad reads | Old processed pool | Same reads re-pooled | Separate pools |
| --- | --- | --- | --- | --- |
| 1 | 3215472 | 11234 | 12837 | 13426 |
| 2 | 3969582 | 13483 | 14898 | 15294 |
| 3 | 4223139 | 13888 | 15246 | 15580 |
| 4 | 3868619 | 13178 | 14915 | 15314 |
| 5 | 4603514 | 14001 | 15339 | 15650 |
| 6 | 4612348 | 13990 | 15278 | 15533 |
| 7 | 4457274 | 13823 | 15204 | 15555 |
| 8 | 3814517 | 12619 | 14445 | 14939 |
| 9 | 3674724 | 12739 | 14462 | 14966 |
| 10 | 3423665 | 12187 | 14102 | 14676 |
| 11 | 3036872 | 11407 | 13340 | 14014 |
| 12 | 1681943 | 6739 | 8734 | 9861 |

The matched ten settings are k15/c1/r1.8,2.55; k17/18/19/c1/r2.55,3.5; and k20/c1/r3.5,5. New counting keeps actual read lengths with ANTK_READ_LIMIT=200. At limit=60, all eleven policies reproduce the earlier output bytes. Global index domain 196596 and the 1050-path cap are retained; block indices are never renumbered locally.

Re-pooling and splitting use exactly the same assigned sequences/qualities, independently checked by a full-record multiset audit. In block12, matched-grid coverage rises from 8,734 re-pooled to 9,861 separated (+1,127), isolating separation itself. The older processed-pool count of 6,739 also changes filtering and length handling. All 319 core settings complete. Five all-raw unfiltered mixed auxiliary settings complete; further work stops after one k18 search takes 700.778 s and overflows 196,322 indices. Incomplete settings are neither zero outcomes nor part of the matched ten.

Across the original 179 settings, separated broad block12 reaches 10,033 versus 6,927 in the old mixed pool: 3,111 gained and five lost, net +3,106. Strict and inclusive matched-ten unions are 9,544 and 9,860, neither above broad's 9,861. At this stage the lowest-budget actual block12 output is k17/c1/r3.5/corrected_vote2: 9,366 correct, 43 wrong, 6,974 erased, worst budget 7,060.

Actual per-block outputs below minimize the measured worst-lane budget. Parameter selection is truth-informed exploratory evaluation; decisions within each file do not read truth. Block12 is updated to the later small-k result.

| Block | k/c/r/policy | Correct | Wrong | Erased | Max budget |
| --- | --- | --- | --- | --- | --- |
| 1 | 17/1/3.5/strict | 13150 | 5 | 3228 | 3238 |
| 2 | 17/1/3.5/strict | 14925 | 5 | 1453 | 1463 |
| 3 | 17/1/2.55/strict | 15170 | 2 | 1211 | 1215 |
| 4 | 17/1/3.5/strict | 15008 | 2 | 1373 | 1377 |
| 5 | 17/1/2.55/corrected_vote5 | 15209 | 4 | 1170 | 1178 |
| 6 | 18/1/2.55/strict | 14988 | 3 | 1392 | 1398 |
| 7 | 18/1/3.5/corrected_vote3 | 15057 | 17 | 1309 | 1343 |
| 8 | 17/1/3.5/strict | 14601 | 9 | 1773 | 1791 |
| 9 | 17/1/3.5/strict | 14611 | 9 | 1763 | 1781 |
| 10 | 17/1/3.5/corrected_vote5 | 14396 | 3 | 1984 | 1990 |
| 11 | 17/1/3.5/strict | 13777 | 4 | 2602 | 2610 |
| 12 | 16/1/3.5/strict | 9973 | 15 | 6395 | 6425 |

## 10. Block12: 174 additional small-k and ratio settings

On the same 1,681,943 broad reads, 174 additional settings bring the block12 total to 353. An initial 129 add k14/k16 and extend k13/k15/k17; twenty refine k16/k15 ratios, and twenty-five refine k14 cutoffs. Input, path limit, index domain and RS policies remain unchanged apart from k/c/r. block12_small_k/plan.json and block12_small_k/parameters.tsv preserve the full grid, including poor outcomes.

| k | Settings | Union at k | Best single candidate | New vs prior 179 |
| --- | --- | --- | --- | --- |
| 13 | 23 | 1603 | 746 | 0 |
| 14 | 67 | 7497 | 5394 | 46 |
| 15 | 59 | 10052 | 8769 | 316 |
| 16 | 50 | 10395 | 10237 | 521 |
| 17 | 39 | 9571 | 9567 | 10 |

The best single candidate count is k16/c1/r3.5: 10,237, including 10,029 with a zero-correction path and 208 requiring inner correction at that setting; r3.75 ties at 10,237. The 353-setting union is 10,659/16,383 (65.0613%), gaining 626 over 10,033. Across all settings, 205 indices require RS correction. There are still 5,724 missing indices; a truth-selected correct-only set is 318 short of the 10,977 erasure-only requirement.

| Actual policy (k16/c1/r3.5) | Correct | Wrong | Erased | Worst 2E+S |
| --- | --- | --- | --- | --- |
| pooled_weight16 | 10077 | 1540 | 4766 | 7826 |
| strict | 9973 | 15 | 6395 | 6425 |

Strict lane errors are [13,14,15,15,15,15], with budgets [6421,6423,6425,6425,6425,6425]. Very short k13/k14 create more collisions and overflow; k16 being optimal does not mean smaller is always better. Re-running the lowest-budget setting without an evaluation reference yields identical bytes for all eleven policies and status/decisions/inner_stats. Earlier read routing remains reference-assisted.

## 11. Full-file decoding after explicit design supplementation

Both simulations retain the actual outputs for blocks 1-11 and inject correct design strands only at missing block12 indices, in ascending index order, never replacing existing strands. Case A retains all 15 wrong strands and adds 1,019, yielding 10,992 correct +15 wrong. Case B first truth-selects the 10,659 successful candidate indices and adds 318, yielding 10,977 correct and no wrong strands. Both candidate selection and injection in case B use truth.

| Case | Injected | Correct | Wrong | Erased | Max budget | ZIP exact |
| --- | --- | --- | --- | --- | --- | --- |
| actual_plus_1019 | 1019 | 10992 | 15 | 5376 | 5406 | True |
| oracle_union_plus_318 | 318 | 10977 | 0 | 5406 | 5406 | True |

The original outer decoder actually processes all twelve blocks in both simulations, producing a 1,268,026-byte ZIP identical to the complete-design baseline and passing every entry CRC. SHA-256 is a7d53ade64bab3f60252e75b28e1741d85fdd7b5b393ffddc95fd141466d8cbb. Injection ledgers and original wrong strands are retained. One fewer injected strand gives budget 5,407; that out-of-bound full decode is not run, only its budget is calculated.

The conclusion is conditional: with these explicitly injected correct strands, the outer implementation recovers the full file. Supplemented simulations and the complete-design baseline are not sequencing-only recovery results.

## Appendix A. Exact grid and execution records

The table lists measured c values for each k/r present in final pooled metrics. Unlisted combinations are not implied to have been tested. Evidence contains individual correctness, errors, overflow, ambiguity, block budgets, commands and timings. Earlier repeats may be superseded by later runs; ledgers retain history, while final metrics count each setting once.

| k | r | Measured c values |
| --- | --- | --- |
| 13 | 1.2 | 9, 49, 99 |
| 13 | 1.5 | 9, 49, 99 |
| 13 | 1.8 | 9, 49, 99 |
| 13 | 2.55 | 1, 3, 9, 19, 49 |
| 15 | 1.2 | 1, 3, 9, 19 |
| 15 | 1.5 | 1, 3, 9, 19 |
| 15 | 1.8 | 1, 3, 9, 19 |
| 15 | 2.1 | 1, 3, 9, 19 |
| 15 | 2.55 | 1, 3, 9, 19, 49 |
| 15 | 3.5 | 1, 3, 9, 19 |
| 17 | 1.2 | 1, 3, 9, 19 |
| 17 | 1.5 | 1, 3, 9, 19 |
| 17 | 1.8 | 1, 3, 9, 19 |
| 17 | 2.1 | 1, 3, 9, 19 |
| 17 | 2.55 | 1, 3, 9, 19, 49 |
| 17 | 3.5 | 1, 3, 9, 19 |
| 18 | 1.8 | 1, 2, 3, 9 |
| 18 | 2.1 | 1, 2, 3, 9 |
| 18 | 2.55 | 1, 2, 3, 9 |
| 18 | 3.0 | 1, 3, 9 |
| 18 | 3.5 | 1, 2, 3, 9 |
| 18 | 5.0 | 1, 3, 9 |
| 18 | 10.0 | 1, 3, 9 |
| 19 | 1.2 | 1, 3, 9, 19 |
| 19 | 1.5 | 1, 3, 9, 19 |
| 19 | 1.8 | 1, 3, 9, 19 |
| 19 | 2.1 | 1, 3, 9, 19 |
| 19 | 2.55 | 1, 3, 9, 19, 49 |
| 19 | 3.0 | 1, 3, 9 |
| 19 | 3.5 | 1, 3, 9, 19 |
| 19 | 5.0 | 1, 3, 9 |
| 19 | 10.0 | 1, 3, 9 |
| 20 | 1.8 | 1, 3, 9 |
| 20 | 2.1 | 1, 3, 9 |
| 20 | 2.55 | 1, 3, 9 |
| 20 | 3.0 | 1, 3, 9 |
| 20 | 3.5 | 1, 3, 9 |
| 20 | 5.0 | 1, 3, 9 |
| 20 | 10.0 | 1, 3, 9 |
| 21 | 1.2 | 1, 3, 9, 19 |
| 21 | 1.5 | 1, 3, 9, 19 |
| 21 | 1.8 | 1, 3, 9, 19 |
| 21 | 2.1 | 1, 3, 9, 19 |
| 21 | 2.55 | 1, 3, 9, 19, 49 |
| 21 | 3.5 | 1, 3, 9, 19 |
| 25 | 2.55 | 1, 3, 9, 19, 49 |
| 31 | 2.55 | 1, 3, 9, 19, 49 |

## Appendix B. Reproduction entry points and script roles

Standalone code archive: File3_code.zip. Preserve its layout and run from its root. It contains the recorded DBGPS2 core, adapters, portable RS-header snapshot, MUSCLE source archive, pinned-version metadata and small validation fixtures. Raw reads are excluded; the explicit fetch stage downloads the matching complete FASTQ and verifies MD5.

```sh
python3 tools/reproduce.py build
python3 tools/reproduce.py smoke
python3 tools/reproduce.py fetch
python3 tools/reproduce.py prepare
python3 tools/reproduce.py graph
python3 tools/reproduce.py msa
```

```sh
python3 tools/reproduce.py inner
python3 tools/reproduce.py partition
python3 tools/reproduce.py small-k
python3 tools/reproduce.py supplement
```

File3 inner analysis requires the original 179 status files; partition needs minimap2 and the inner build; small-k requires exported broad block reads; supplementation needs the small-k report and actual first-eleven-block outputs. Package README documents dependencies. Union reports need NumPy, plots additionally need Matplotlib; Python 3.12 is recommended. Full replay includes hundreds of millions of reads and substantial runtime.

| Script | Role |
| --- | --- |
| download.py / fetch_sources.py | Verified downloads and source pins |
| prepare_reads.c / prepare_references.py | Length preprocessing and reference transcoding |
| graph_adapter.c / run_graphs.py | Strict traversal and parameter grid |
| msa_baseline.py / inner_select.cpp | MSA and common inner-code decisions |
| rs_driver.cpp / fast_gf14.hpp / legacy_random.hpp | Outer decode, field acceleration and legacy sampler |
| graph_inner_rs.cpp / run_inner_graph.py | Candidate correction and eleven policies |
| run_block12_analyzer.py / report_block12_analyzer.py | Independent coverage and graph diagnostics |
| split_reads_by_alignment.py / split_aligned_blocks.cpp | Full alignment and sequence-preserving routing |
| run_block_decode.py / report_block_decode.py | Same-read re-pooling and block comparison |
| run_block12_small_k.py / report_block12_small_k.py | 174 added settings and 353-setting union |
| simulate_file3_supplement.py | Explicit truth injection and full-file verification |

## Appendix C. Limitations, sources and evidence map

All best settings are exploratory choices on this dataset, with no held-out validation or multi-seed confidence intervals. MSA is not a complete LSH reproduction. Coverage, candidate unions, actual decisions and file recovery are distinct outcomes. Filtering, length handling, partitioning and added truth information are disclosed separately. Unattempted, interrupted and diagnostic failures are not represented as completed measurements.

Primary local evidence is under evidence/File3/: original summary.json, graph_metrics.json, graph_parameters.csv, run_ledger.json, verified-input records and decode validation. MANIFEST.sha256 records delivered-file digests; EVIDENCE_MAP.json links each evidence item to its original workspace path.

[Article: DOI 10.1038/s41467-020-19148-3](https://www.nature.com/articles/s41467-020-19148-3). [Primary repository](https://github.com/MLI-lab/noisy_dna_data_storage), pinned to 41422be448b3c1da9041c3a2967d2754d2ad8b1a. [RS repository](https://github.com/reinhardh/dna_rs_coding), checkout 455e1a5182bf64281f54d434a4ed5ca456c2373e, historical headers from 8e7e3b5. [Later designs](https://github.com/fml-ethz/dt4dds-challenges_notebooks), pinned to 861d3b95b52d1b121e9d85040f20c8ecf0cda2d4.
