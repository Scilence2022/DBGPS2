# Fragment Recovery and Complete Decoding Report: DBGPS2 on the DNA Fountain Dataset (English Version)

- **Project**: DBGPS2 (de Bruijn graph path-search assembler) × DNA Fountain dataset
- **Paper**: Erlich & Zielinski, *DNA Fountain enables a robust and efficient storage architecture for digital data*, Science 355(6328):950–954, 2017 (doi:10.1126/science.aaj2038)
- **Sequencing data**: ENA ERR1816980 (BioProject PRJEB19305, master pool, Illumina MiSeq PE151, 16,692,317 read pairs = 33.4M reads; R1 md5 `f53b849e…`, R2 md5 `acd4269a…`, both matching the official ENA checksums)
- **Paper code**: github.com/TeamErlich/dna-fountain (Python 2)
- **Workspace**: git worktree `DBGPS2-dnafountain`, branch `dnafountain-adapter`
- **Bottom line**: from the complete data DBGPS2 assembled **71,990 of the 72,000** designed strands (99.99%); using only its output, the fountain code was fully decoded and the recovered file's **md5 = `8651e90d3a013178b816b63fdbb94b9b` matches the paper's original input bit-for-bit (error-free recovery of 2,146,816 bytes)**

> **Correction (2026-11-01).** Earlier versions of the paper-style comparator had two faults: its RS filter accepted words in which only the parity needed correction, and its biochemical screen tested 136 nt although the encoder screens the full 152-nt block, which wrongly removed about 5,000 designed strands. The comparator was fixed and rerun on the raw reads, and the numbers here are updated. The earlier claim that DBGPS2 recovers 5,007 strands the read-level pipeline cannot is withdrawn.

---

## 1. Analysis of the Paper's Coding Scheme

### 1.1 Encoding pipeline

1. The payload (a 2,146,816-byte compressed tarball) is split into **67,088 segments of 32 bytes**.
2. **Luby transform** (robust soliton degree distribution, c=0.025, δ=0.001): each droplet is the bitwise XOR of a random degree-d subset of segments.
3. **Droplet layout** (38 bytes = 152 nt): `4B seed ‖ 32B payload ‖ 2B RS(38,36) parity`
   - Seed = sequential states of a 32-bit Galois LFSR (polynomial 1+x²⁵+x²⁶+x³⁰+x³², initial state 0b001010101);
   - RS = reedsolo `RSCodec(2)`: GF(256), primitive polynomial 0x11D, α=2, fcr=0, systematic (data first, parity last).
4. **Bytes→DNA**: `{00,01,10,11}→{A,C,G,T}`, MSB first — **identical to DBGPS2's built-in mapping**.
5. **Biochemical screening**: homopolymer runs ≤3, GC 45–55% (applied to the full 152-nt sequence including RS); the LFSR advances and rejected droplets are retried. Empirically the design consumed roughly the first **577,000 LFSR states** to emit **72,000 oligos** (~12.5% acceptance).
6. Final oligo = 26-nt 5′ adapter (`GTTCAGAGTTCTACAGTCCGACGATC`) + 152-nt data block + 21-nt 3′ adapter (`TGGAATTCTCGGGTGCCAAGG`).

### 1.2 The paper's original decoding pipeline (our comparator)

PEAR-merge read pairs → keep 152-nt inserts → abundance sort (sort | uniq -c) → feed droplets one by one to the decoder, which applies **RS detect-only** semantics (`--max_hamming 0`), a 136-nt biochemical screen and seed de-duplication, then runs a peeling decoder to invert the Luby transform. Reported: with full data the decoder consumed **69,870 valid oligos** to finish recovery; at 750k reads (one tile) the dropout was 1.3% with 20/20 perfect decodes.

### 1.3 Key empirical facts about the ENA data

- The submitted FASTQ files are **already adapter-trimmed**: R1 = block[0:151] and rc(R2) = block[1:152] (verified on 20,000 sampled pairs: 52.8% satisfy the relation exactly; the rest are PhiX / N-containing / erroneous).
- Library reads make up ~47–53% of the files; the rest is PhiX spike-in.
- Both mates contain 16,692,317 records each when intact (strictly paired).

---

## 2. DBGPS2 Adaptations (implemented on the worktree branch)

| Difference | Handling |
|---|---|
| Base mapping, 16-nt index, fixed 128-nt payload, RS(2B) frame | Natively compatible (`-i 16 -l 128 -E rs:2B`); the C RS implementation produces byte-identical parity to reedsolo (cross-validated, including a 1-byte correction round trip) |
| Reads begin at the seed; no left flank for anchoring | Option A: k=i=16 with no primer anchor; **final choice**: prepend a constant 10-nt prefix to every R1 read + `-A CCTGCAGAGT`, k=21 |
| Sparse 32-bit LFSR seeds (not a contiguous -a/-b range) | New `-I <index-list file>` option (list extracted from R1 against the legal LFSR state set) |
| Paper semantics max_hamming=0 | New environment variable `DBGPS_DETECT_ONLY=1`: RS detect-only, no correction |
| End-of-search 3′ one-base lookahead (see §4.6) | Prepend 'A' to R2 reads (the reverse-complement side then carries 'T', the first 3′-adapter base) |
| Fixed coverage ratio vs Poisson coverage | Parameter tuning: `-c 0 -r 6.0` optimal |

After adaptation, `make test` passes all regression checks (default behaviour unchanged).

---

## 3. Data Acquisition (including failures and repairs)

1. ENA throttles a single connection to ~12–15 KB/s; per the user's instruction all downloads were **single-threaded**. R1 (1.476 GB) completed after ~6 h with an md5 matching the official one.
2. **R2 was corrupted on first download**: interrupted resumes glued a mid-file region that stays deflate-valid but decodes to garbage (intermittent zones; first bad record at ~7,884,712; file md5 mismatch, gzip trailer CRC error). Bisection showed only 189 decompressed bytes missing at the tail, but in-place byte damage mid-file. At that time 7,828,857 clean pairs (~47% of the run) were salvaged to complete the experiments.
3. **R2 re-download** (at the user's request): single stream at ~441 KB/s, ~51 minutes; **md5 = `acd4269a5e581a9d1c3117dfa6ad6fc7` exactly matches ENA**, `gzip -t` passes; replaced the corrupted file. The intact R2 contains 16,692,317 records (the corrupted copy was missing 764,300).
4. Final experiments are based on **both mates complete and md5-verified**.

---

## 4. Full Log of Decoding Attempts (chronological, including failures and root causes)

### 4.1 Smoke test (partial data, 145k reads)

- Seed extraction: first 16 nt of R1 filtered against a 600k-state LFSR window → 45,800 valid seed observations (random PhiX hits rejected).
- DBGPS2 `-E none` (no checksum, unique-path requirement): 5,512 strands of 144 nt; 79.6% exactly matched the naive read-level consensus.
- Validation trick: `-E none` output lacks the RS bases → derive the 8-nt parity forward from the 36 data bytes with reedsolo and compare with the read majority vote at block[144:151] (a 14-bit independent check); 4,205/5,512 passed.

### 4.2 LFSR window bug discovered

An initial 400k-state filter showed observed seeds uniformly spread over 0–400k with the boundary suspiciously close to 400k. Extending to 2M states revealed the true boundary ≈ 577k (the number of states consumed by the design). Fixed with a 600k window (~1,056 expected random hits, harmless).

### 4.3 First rs:2B run and the "35% correct" mystery (375k pairs)

With `-c 0 -r 2.55` DBGPS2 produced 13,402 strands, but only 35% matched the consensus, with Hamming distances up to 40–90. Investigation, layer by layer:

1. **Our own verifier bug**: verify_strands compared the 36-byte corrected data against the full 38-byte word, so the "corrected" counter was always true (cosmetic).
2. **Cross-oligo chimeras**: at k=16 the pool contains ~7,000 shared 16-mers; the DFS can jump from oligo A onto real sequence of oligo B and assemble **another valid codeword** (indistinguishable by RS; the 1,589 screen-violating outputs were exactly of this kind). → k=21 (shared-kmer probability vanishes) eliminated them.
3. **Naive majority consensus is unreliable**: testing whether the consensus itself forms a codeword showed 90% of mismatched consensuses do **not** — mixed PCR/synthesis chimeric subpopulations corrupt the per-position vote, while the codeword subpopulation selected by DBGPS2+RS is the true design. **RS codeword validity is the stronger criterion.**
4. **RS(2) miscorrection**: candidates with 2 byte errors get "corrected" onto nearby wrong codewords → added `DBGPS_DETECT_ONLY=1` mirroring the paper.

### 4.4 Coverage-ratio tuning (750k pairs)

The fixed `-r 2.55` leaves a ~0.02% 136-step survival at λ≈6 (Poisson noise hits the window); the successes were high-coverage hotspots. Sweep results:

| Parameters | Yield |
|---|---|
| k=16, r=2.55 (initial) | 13,402 (35% correct) |
| k=31 (prefix anchor), r=2.55 | 8,843 (49.8% correct) |
| k=21, r=2.55 | ~13k |
| **k=21, -c 0, r=6, detect-only** | **27,626 (99.66% agreement with the paper-style pipeline)** |

One **run-to-run nondeterminism** was also observed (same command: 12,618 vs 23,809; four reruns including single-threaded were stable at 23,809) — not reproducible, recorded as an open question for upstream.

### 4.5 Fountain decoder port and debugging (Python 3 re-implementation)

- Strict Python 2 semantics: `random.seed(int)` uses init_by_array in both versions; `random.sample` implements the **float-multiplication path of the actual 2.7.18 source** (set branch `j=int(random()*n)` with rejection, unsorted; pool branch likewise); degree sampling and the robust-soliton CDF compared expression by expression. The 2.7.18 source itself was fetched from python.org for line-by-line verification.
- **First decode stalled (194/67,088)**: initially suspected a PRNG porting error; a synthetic round-trip once yielded 16/2,000 — finally traced to a flaw in my own synthetic test (leading zero bytes of small integer seeds → runs of 12 A's violating the homopolymer screen; the first pass occurred exactly at tries=16,777,217=0x01000000, a smoking gun).
- **The real bug**: the peeling decoder never registered droplets into chunk_to_droplets (glass.pyx's addDroplet registers before updateEntry), so cascades could not happen. Fixed: 194 → 1,539 solved.
- **Density threshold**: synthetic curves at K=67,088 show 94.2% density solves only 2.3% (exactly matching the real data's 1,539/67,088 = 2.3%, which independently confirms the PRNG port); the transition sits at **~104% (69,771 droplets → 100%)**, perfectly consistent with the paper's 69,870.
- Also discovered: the **decoder-side biochemical screen (136 nt) differs from the design-side screen (152 nt)** — the design only guarantees the latter; ~7% of true droplets fail the decoder-side screen. For consensus droplets already validated as full 38-byte RS codewords, decoding uses `--no_screen` (its noise-filtering purpose is moot for consensus output).

### 4.6 "Adding R2 destroyed 46k strands" and the 3′ one-base lookahead (full data)

After cleaning (7.83M pairs): the R1-prefixed file alone assembled 68,275 strands; adding R2 collapsed the yield to 21,767. Investigation:

1. A C probe read the hash table directly: every step of the true path lies inside [parent/6, parent×6] — ratio pruning excluded.
2. Frontier tracing: the frontier stayed at 1 through step 136 and **died at step 137**. Root cause: `while(p <= path_len)` makes the search consume **one base beyond** the nominal path — the 3′ one-base lookahead (the code comment itself says "the terminal node is the historical one-base primer lookahead"). The original DBGPS2 data keeps 3′ primers in reads, so the lookahead always exists; adapter-trimmed Fountain reads lack it → guaranteed death (only randomly covered exceptions survive).
3. **Fix**: prepend 'A' to R2 reads (rc side gains 'T', the first 3′-adapter base). After verifying on a single seed: **21,767 → 71,983 strands**.

### 4.7 End-to-end success

- 71,983 droplets (47% of data) decoded with `--no_screen`: 67,088/67,088 chunks solved, **md5 matched exactly on the first complete attempt**.
- After the R2 re-download, the definitive run on complete data is reported in §5.

---

## 5. Final Results (complete 33.4M reads, both mates md5-verified)

| Metric | Paper | Paper-style pipeline (our reproduction) | **DBGPS2** |
|---|---|---|---|
| Designed / encoded strands | 72,000 | — | — |
| Valid droplets, full data | decoder consumed 69,870 | 71,998 (all exact) | **71,990 (99.99%, all exact)** |
| 750k pairs (one-tile equivalent) | dropout 1.3%, 20/20 perfect | 70,283 (decodes, md5 matched) | 63,327 (94.4%, below the 104% decoding threshold; 1,466/67,088 chunks) |
| Cross-agreement with paper-style | — | — | 71,988 common, **all identical**; union covers all 72,000 designs |
| Exclusive recoveries | — | 10 | 2 |
| End-to-end recovery | perfect | 67,088/67,088 chunks, md5 bit-exact | **67,088/67,088 chunks, md5 bit-exact** |
| Runtime | decoding ~9 min (2016 laptop) | ~50 min (pure Python) | k-mer counting ~2 min + **assembly 4.8 s** + decoding ~1 min |

### 5.1 Direct answers to the three questions

1. **Strands encoded by the paper**: 72,000.
2. **Strands decoded by the paper**: 69,870 (valid oligos consumed by the time the decoder finished).
3. **Strands assembled by DBGPS2**: 71,990 (all correct, certified by the end-to-end md5).

### 5.2 Conclusions

- **Applicability**: the DNA Fountain coding scheme (mapping / index / fixed payload / RS frame) is **fundamentally compatible** with DBGPS2's structural assumptions; only two small adaptations (`-A`, `-I`), the k=21 prefix anchor, detect-only semantics and the R2 'A'-prepend were needed.
- **Capability**: DBGPS2's per-strand accuracy is ~100% (double-verified by RS codeword checks and the end-to-end md5); one assembly pass on full data covers all redundancy needed for decoding. On this high-quality library the read-level pipeline does at least as well at full depth (71,998 against 71,990 strands).
- **Boundary**: at very low coverage (one-tile equivalent, λ≈11) "waiting for one perfect read pair" (paper-style) clearly beats "graph consensus + ratio pruning" (DBGPS2); the fixed ratio window is conservative under Poisson coverage.

### 5.3 Improvement suggestions for the DBGPS2 project

1. Merge `-I` (sparse index list) and `-A` (custom primer) into mainline — external datasets generally need them.
2. Provide an off switch (or at least documentation) for the 3′ one-base lookahead — on adapter-trimmed datasets it silently kills every search, which is very hard to diagnose.
3. Consider noise-adaptive coverage-ratio windows (fixed windows inevitably prune true paths at low Poisson coverage).
4. Investigate the one unreproduced run-to-run nondeterminism.
5. RS correction with only 2 parity bytes carries miscorrection risk; a detect-only option is worthwhile (implemented on this branch).

---

## 6. Reproduction

```bash
# 0) Data: ENA ERR1816980 (ERR1816980_{1,2}.fastq.gz), verify md5
# 1) Convert to FASTA (plain R1 / prefixed R1 / plain R2 / 'A'-prepended R2)
gzip -dc ERR1816980_1.fastq.gz | awk 'NR%4==1{print ">" NR} NR%4==2{print}' > r1.fa
gzip -dc ERR1816980_1.fastq.gz | awk 'NR%4==1{print ">" NR} NR%4==2{print "CCTGCAGAGT" $0}' > r1p.fa
gzip -dc ERR1816980_2.fastq.gz | awk 'NR%4==1{print ">" NR} NR%4==2{print "A" $0}' > r2A.fa
# 2) Seed extraction (LFSR-filtered)
python3 analysis/extract_seeds.py --r1 r1.fa --out seeds_all.tsv
awk '{print $1}' seeds_all.tsv > seedlist.txt
# 3) DBGPS2 fragment recovery
DBGPS_DETECT_ONLY=1 ./DBGPS2 -k 21 -i 16 -l 128 -E rs:2B -A CCTGCAGAGT \
    -c 0 -r 6.0 -I seedlist.txt -o out r1p.fa r2A.fa
# 4) Fountain decoding and verification
cut -f2 out.strands > droplets.txt
python3 analysis/decode_fountain.py droplets.txt --no_screen --out decoded.bin
md5 decoded.bin   # expected 8651e90d3a013178b816b63fdbb94b9b
```

Dependencies: `python3` (`reedsolo` ≥1.7), `gzip/awk`, the adapted `DBGPS2` (patch included as `patch/DBGPS2-dnafountain.patch` in the package).

## 7. Artifact locations

- Final assembled strands: `analysis/runs/fullint/dbgps_do.strands` (71,990 × 152 nt)
- Final decoded file: `analysis/runs/fullint/decoded_do.bin` (md5 8651e90d…)
- Comparator pipeline output: `analysis/runs/fullint/paper_droplets.txt`
- Per-stage run directories: `analysis/runs/{onetile,tile750k,full,fullint}/`
- Code package: `DBGPS2-DNAFountain-pipeline.tar.gz`
