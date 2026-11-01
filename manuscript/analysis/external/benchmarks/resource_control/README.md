# Single-worker and three-worker reconstruction controls

This new experiment retains the complete DNA Fountain, Organick File 16 and HT176 inputs, all nine configurations and fixed search/Bloom parameters of the preceding eight-worker matrix. Each of the two resource settings has three complete computational repeats: 162 new measurements. Previous timing values are not reused.

## Resource protocols

- **Serial:** one graph-pipeline worker performs read, extract and insert stages sequentially; `kt_for(1)` runs its counting callback inline on that same worker. Seed enumeration, query-filter construction and path search are serial. The MSA program directly maps over eligible clusters without a Python process pool; it waits for each single-threaded MUSCLE3 child, then votes before starting the next cluster. This controls algorithmic parallelism to one active reconstruction worker; it does not bind that worker to a particular physical CPU core.
- **Three primary workers:** three count insertion workers plus the existing three-stage pipeline; three Python MSA workers, each waiting for one single-threaded MUSCLE3 child and performing its voting between calls. Both methods retain a coordinator and other pipeline/IPC overhead. This matches the requested primary worker counts and is **not a strict three-core or three-total-thread budget**.

Measured reconstruction jobs run sequentially on the same workstation, with cyclic method order and alternating resource-setting order. Numerical-library thread limits are one. Background user applications are not controlled; no affinity, dedicated idle host or cold-cache regime is claimed. Pre-run process snapshots are recorded. The OS process count is not treated as a measure of active CPU parallelism.

## Unchanged reconstruction rules

The original saved counter, Bloom APIs, path acceptance, read selection, MUSCLE3 version, voting, format validation and conflict rejection are retained. Benchmark-only source changes expose count/pipeline sizes and add getrusage stage instrumentation. Serial MSA bypasses process-pool dispatch; three-worker MSA retains the original process-pool implementation. No production files were edited or committed.

Inputs and parameters are inherited from ../speed_extension. The expanded matrix is no Bloom at two c values, query B8/B16, counting F32/F34 and combined B12/F34, plus the adapted DBGPS absolute-cutoff rule and exact-prefix/MSA/vote comparator. The DBGPS comparator emulates the rule in the same saved core and is not the original 2022 binary. HT176 is the local 176-day sample; the specific sample's public accession remains unverified. No HT5 speed runs are included.

## Timing

End-to-end wall time starts before subprocess launch and ends after native reconstruction, validity filtering, reference-free conflict rejection and sorted output writing. Ground-truth scoring, output-hash comparison, outer decoding and preprocessing occur outside the timer. Completed outputs must match the corresponding previously verified eight-worker configuration; that is a post-timing correctness assertion, not reference-based reconstruction.

CPU time is user plus system time from `/usr/bin/time -l`, plus calling-harness self CPU over the corresponding interval (including output canonicalization). Native user/system and harness user/system are retained separately. Graph phase CPU uses `getrusage(RUSAGE_SELF)` over count, filter construction, seed and search, covering its native threads. MSA phase CPU combines self and reaped-child usage after each phase; the pool is closed and joined before phase-end sampling. Clustering, alignment/voting and validation are separate stages. Stage wall-time sums omit process startup/shutdown, graph destruction and final harness output work; these remain included in complete reconstruction. Native peak RSS is reported for graph executables; MSA process-tree RAM is not inferred from per-process maximum RSS.

`source_manifest.json` records saved source/executable hashes. `smoke_validation.json` records 50 passed format, scheduling, Bloom and fragmented-read checks. `resource_results.json` updates only after complete successful measurements and resumes from complete per-run result records. All raw stdout/stderr and stage records are retained.

The portable package also retains the pinned MUSCLE 3.8.31 source archive and the original native build script/provenance. Input/reference paths must be configured for another host; large sequencing FASTA/FASTQ files remain outside the archive.
