# Optional Bloom filters

Both assemblers support query filtering and two-pass singleton removal. Neither
is enabled unless requested. This document records the implementation selection
and validation on 2026-10-15.

## GitHub implementation comparison

The integration targets are C, canonical packed integer k-mers, the existing
`khashl` partitions, GCC/Clang, and both x86-64 and ARM64. A small scalar backend
fits these constraints without changing the build or adding external libraries.

| Repository | Relevant properties | Fit for this project |
| --- | --- | --- |
| [lh3/yak](https://github.com/lh3/yak) | MIT; small C implementation; 64-byte blocked bitmap; integer hashes; partitioned two-pass k-mer counting | Selected: matches this project's rolling hash, partition workers, and `khashl` conventions closely |
| [jvirkki/libbloom](https://github.com/jvirkki/libbloom) | BSD-2-Clause; general C API for byte buffers; MurmurHash2; parameterized capacity/error rate | Good general library; adapting an already packed integer key and partitioned counting requires more surrounding code than yak |
| [jbapple/libfilter](https://github.com/jbapple/libfilter) | Apache-2.0; C and C++ APIs; blocked/static/growing filters; SIMD-aware block implementation | Useful broader toolkit, with more code and platform paths than this integration needs |
| [PeterRK/PageBloomFilter](https://github.com/PeterRK/PageBloomFilter) | BSD-3-Clause; page-based filter; native C/C++ fast path restricted to little-endian x86-64 and AArch64 | Plausible optimized alternative; imposes a narrower architecture contract and a larger integration surface |

This is a fit assessment, not an upstream performance ranking. No claim is made
that yak is faster than every alternative. The relevant primary sources are
yak's [bbf.c](https://github.com/lh3/yak/blob/ef033d720e229593b32ad454af14ccea2120ec13/bbf.c),
[htab.c](https://github.com/lh3/yak/blob/ef033d720e229593b32ad454af14ccea2120ec13/htab.c),
and [count.c](https://github.com/lh3/yak/blob/ef033d720e229593b32ad454af14ccea2120ec13/count.c),
libbloom's [API/source](https://github.com/jvirkki/libbloom/blob/master/bloom.c),
libfilter's [block header](https://github.com/jbapple/libfilter/blob/master/c/include/filter/block.h),
and PageBloomFilter's [platform documentation](https://github.com/PeterRK/PageBloomFilter#c).

## Provenance and local changes

`dbgps_bloom.h` adapts yak's `bbf.c`, definitions from `yak.h`, and the full-width
integer mixer from `yak-priv.h`, pinned to commit
`ef033d720e229593b32ad454af14ccea2120ec13`. The upstream MIT notice is retained in
the header. No network fetch, separate installation, C++ runtime, SIMD flags, or
new link dependency is required.

Local changes include checked allocation and size arithmetic, automatic query
filter sizing, portable 64-byte alignment, and a read-only membership operation.
The graph hash contains only `2*k` useful bits, so the backend remixes it across
64 bits before choosing a block and probe positions. This matters for small k.
Four probes touch one aligned 512-bit block. Queries do not modify filter bits.

`dbgps_bloom_graph.h` shares the graph integration between DBGPS2 and DBGPS2-RPC.
The existing partition workers give each first-pass filter a single writer.

## Query option: `-B INT`

The default `-B 0` allocates no query filter. Values 1 through 64 specify a minimum
bitmap budget in bits per distinct retained canonical k-mer. `-B 12` is a useful
starting value. The bitmap rounds up to a power of two, with a 64-byte minimum;
it therefore uses approximately `B*N/8` to `2*B*N/8` bytes, plus small metadata
and alignment overhead. An empty graph still has the minimum bitmap when enabled.

The filter is built from the exact table after counting. A negative probe returns
zero; a positive probe continues to the exact table, which supplies the real
coverage, including zero on a Bloom false positive. No singleton is removed by
`-B`. The same path covers assembly and JSON-RPC `kmer/cov` calls. This adds
memory and build work, and helps most when absent lookups dominate; small graphs
or queries rich in present keys may be slower.

## Counting option: `-F INT`

The default `-F 0` retains ordinary one-pass counting. Values 23 through 40
specify the base-2 logarithm of the total first-pass bitmap bits across all hash
partitions. Examples:

| Option | Total bitmap |
| --- | ---: |
| `-F 23` | 1 MiB |
| `-F 27` | 16 MiB |
| `-F 30` | 128 MiB |

The budget is for filter bits, not total process memory. Metadata, per-partition
alignment, input buffers, and the exact candidate table also consume memory.
An undersized or saturated filter admits more singleton candidates and can
eliminate the memory saving, while the final counts remain exact.

1. Scan every input file, using canonical hashes. First sightings only set bits;
   probable repeats enter the exact candidate table with count zero.
2. Release the first-pass Bloom bitmaps, replay all input files, and increment
   only exact candidate keys, using the existing saturating 14-bit counters.
3. Remove candidates whose exact coverage is below two. This also removes
   singletons admitted by Bloom false positives.
4. Build the optional `-B` query filter from the retained exact graph.

Every k-mer occurring at least twice is selected and counted exactly, including
repeats spanning files, pipeline blocks, and opposite strands. Singletons are
always absent from the final graph. Coverage-one queries become zero, and graph
metadata counts fewer nodes. Paths that require coverage-one nodes can change.
The existing `-c` threshold and ratio logic remain in effect after filtering.

Both scans use the same read-length restriction, ambiguity handling, and mapping
configuration as ordinary counting. All inputs must be readable regular files,
including compressed files, and remain unchanged between scans. FIFOs and
devices are rejected before opening. Malformed FASTQ and truncated gzip inputs
return an error instead of publishing an incomplete filtered graph.

## Examples and API

```bash
# Query filter only: identical graph and coverage semantics.
./DBGPS2 -B 12 -a 1000 -b 1001 reads.fq.gz

# Singleton removal only: lower graph memory when singleton noise dominates.
./DBGPS2-RPC -F 27 -a 1000 -b 1001 reads.fq.gz

# Both, including server mode.
./DBGPS2-RPC -S -B 12 -F 27 reads1.fq.gz reads2.fq.gz

# Explicitly keep both disabled.
./DBGPS2 -B 0 -F 0 -a 1000 -b 1001 reads.fq.gz
```

The original `dbgps_load(...)` API keeps its defaults. For two-pass counting,
use `dbgps_load_filtered(..., count_bloom_shift)`, whose final argument is zero
or the `-F` exponent. `dbgps_set_bloom(g, 12)` builds/replaces the query filter;
zero disables it. It returns -1 on error while preserving the previous filter.
`dbgps_bloom_bytes(g)` reports allocated query bitmap bytes.

API callers must stop queries while configuring filters or loading data. The
server uses a read/write lock: ordinary queries share read access, and
`graph/load` holds exclusive access until table and filter updates finish.
`graph/info` includes `bloom.enabled`, `bits_per_key`, `bitmap_bytes`,
`count_shift`, and `singletons_removed`. Here `enabled` describes the query
filter; `count_shift` independently describes the counting mode.

With `-F`, `graph/load` rebuilds from the recorded input filenames plus the new
files. This restores historical singletons that now occur twice, avoids lost
first occurrences, and preserves the entire old graph if rebuilding fails.
Keep those input files available and unchanged. Reload time grows with input
history, and peak reload memory includes the old and replacement graphs.

Without `-F`, additional files append counts as before, then rebuild an enabled
query filter. A failed append can retain partial counts; the stale filter is
disabled so exact queries remain usable. Its requested bit budget is retained
for rebuilding on the next successful append.

## Validation and synthetic benchmark

`make test` covers both tools and both Bloom modes, exact fallback under forced
false positives, singleton removal, count saturation, Ns, reverse complements,
cross-file and cross-block repeats, multiple thread counts, all codec regression
cases, replay after server loads, concurrent readers during loads, allocation
failure cleanup, malformed FASTQ, truncated gzip, FIFO rejection, and explicit
disable options. The complete suite passed both normally and with
`-fsanitize=undefined -fno-sanitize-recover=all`. A deterministic primitive test observed 132 positives among
100,000 absent keys at a requested 12 bits/key, with no false negatives; this is
a test fixture result, not a promised false-positive rate.

To reproduce the optional benchmark outside the checkout:

```bash
mkdir -p /tmp/dbgps-bloom-benchmark
cc -O3 -Wall -I. -o /tmp/dbgps-bloom-benchmark/bench tests/bench_bloom.c \
  dbgps_core.c dbgps_codec.c kthread.c -lz -lpthread
/tmp/dbgps-bloom-benchmark/bench generate /tmp/dbgps-bloom-benchmark/reads.fa
/usr/bin/time -l /tmp/dbgps-bloom-benchmark/bench /tmp/dbgps-bloom-benchmark/reads.fa 0 0
/usr/bin/time -l /tmp/dbgps-bloom-benchmark/bench /tmp/dbgps-bloom-benchmark/reads.fa 12 0
/usr/bin/time -l /tmp/dbgps-bloom-benchmark/bench /tmp/dbgps-bloom-benchmark/reads.fa 0 23
/usr/bin/time -l /tmp/dbgps-bloom-benchmark/bench /tmp/dbgps-bloom-benchmark/reads.fa 12 23
```

Use `/usr/bin/time -v` on Linux instead of macOS's `-l` flag. This fixture has
1,000,000 distinct canonical 31-mers: 900,000 occur once and 100,000 occur three
times, for 1,200,000 reads. It performs 10,000,000 queries, approximately 95%
absent, with a 100,000-base pipeline block. This deliberately isolates graph
storage from the CLI's larger default input buffers. Reported query time also
includes generation and canonicalization of query keys.

A representative run on this macOS 26.4.1 ARM64 host, with Clang `-O3`:

| Query / counting options | Retained k-mers | Build seconds | Query ns/op | Peak RSS MiB |
| --- | ---: | ---: | ---: | ---: |
| Both disabled | 1,000,000 | 0.203 | 55.70 | 35.00 |
| `-B 12` | 1,000,000 | 0.212 | 29.81 | 36.78 |
| `-F 23` | 100,000 | 0.403 | 46.96 | 11.34 |
| `-B 12 -F 23` | 100,000 | 0.390 | 26.86 | 11.42 |

All four query checksums were 1,500,000. The query-only option cut query time by
about 46%; counting-only cut peak RSS by about 68% while approximately doubling
build time. These single-host synthetic results do not establish a speedup or
memory reduction for real sequencing inputs, other singleton rates, hit rates,
pipeline sizes, or hardware.

AddressSanitizer could not run on this host: even an empty ASan executable hung
inside the system allocator/dyld initialization, before `main`. The documented
`make test asan=1` target remains available for hosts with a working runtime.
