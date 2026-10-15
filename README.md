# DBGPS2

DBGPS2 is a collection of tools for DNA data storage analysis and processing, with a focus on sequence assembly, cross-link detection, and coverage analysis.

## Overview

This project provides software tools for analyzing and reconstructing DNA sequences from data storage systems. The main components are:

- **DBGPS2**: A read-length constrained sequence assembler with dynamic coverage ratio support (formerly DBGPS-dy2-rl)
- **DBGPS3**: Next-generation de Bruijn graph assembler with interactive server mode support

## Installation

### Prerequisites

- GCC compiler
- zlib development libraries

### Building the software

Clone this repository and build using make:

```
git clone https://github.com/your-username/DBGPS2.git
cd DBGPS2
make
```

To build individual components:

```
make DBGPS2
make DBGPS3
make dbgps3-client
```

To build all targets:

```
make all
```

### Regression checks

Run `make test` (requires Python 3) to check missing k-mer lookups, bounded path
searches, CRC ambiguity, allocation cleanup, and exact strand reconstruction
across repeated CLI runs, thread counts, coverage cutoffs, and ratio sweeps.
They also check all 72 supported codec block widths, custom mapping reconstruction,
forbidden words, CRC padding, malformed codebooks, and larger indices. DBGPS3's
live JSON-RPC checks cover mapping metadata, exact large indices, invalid layouts,
and the interactive client's index shortcut.
Test binaries and generated reads are stored in temporary directories.
Use `make test asan=1` to run the same checks with AddressSanitizer.

## Tool Descriptions

### DBGPS2

A read-length constrained sequence assembler with dynamic coverage ratio handling.

#### Features
- Skip option for skipping initial ratio checks at specified positions
- Output prefix option for log and strands files
- Ratio range iteration support
- Improved path finding algorithm
- Configurable fixed-block DNA/bit mapping with byte lookup acceleration

#### Usage
```
DBGPS2 [options] <input file>
                  [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]

Options:
  -k INT     k-mer size
  -i INT     length of index
  -l INT     data encoding length
  -t INT     number of threads
  -c INT     k-mer coverage cut-off for exclusion of noise k-mers
  -C INT     Switch on k-mer coverage testing mode
  -a INT     Initial index
  -b INT     End index
  -r FLOAT   Maximum coverage ratio
  -R FLOAT   Upper bound for ratio range iteration (0 = disabled)
  -s INT     Number of initial positions to skip ratio checks for
  -o STR     Output prefix for log file and strands file
  -M FILE    Custom DNA/bit mapping (default: none; built-in A=00,C=01,G=10,T=11)
```

#### Custom base-to-bit mappings

`DBGPS2 -M codebook.map` supports a reversible, context-independent mapping from
`N` bases to `M` bits, with `1 <= N <= 8` and `1 <= M <= 16`, `M <= 2N`.
Without `-M`, DBGPS2 uses the same built-in encoding as its k-mers:
`A=00, C=01, G=10, T=11`. Index bits are used directly as packed DNA, and CRC
validation reads the reconstructed bytes directly, without a mapping table or
an extra conversion. For example, the byte `00011011` is the DNA word `ACGT`.
DBGPS2 and DBGPS3 share this default and the same codebook format.

**Default behavior change:** strands encoded with the original special two-base
permutation now require `-M examples/legacy-2base.map` explicitly.

The text format specifies the block widths followed by `DNA BITS` rows. For example:

```text
# One base to two bits; bit order is most significant first.
bases 1
bits 2
A 00
C 01
G 10
T 11
```

Both headers must precede the rows. Blank lines, `#` comments, arbitrary row
order, and lowercase DNA are accepted. All `2^M` bit words must appear exactly
once, each assigned to a different DNA word. Unlisted DNA words are forbidden;
candidates containing them fail validation. Duplicate, incomplete, or malformed
codebooks fail before k-mer counting. Variable-length codes, multiple DNA words
for the same bit word, and mappings that depend on preceding bases are outside
this format.

Ready-to-use examples:

- [legacy-2base.map](examples/legacy-2base.map): explicitly enables the original two-base permutation.
- [acgt-1base.map](examples/acgt-1base.map): an explicit codebook equivalent to the default built-in encoding.
- [3base-5bit.map](examples/3base-5bit.map): an example of 32 allowed DNA triplets.

```bash
make DBGPS2
./DBGPS2 -a 1000 -b 1001 -o assembled reads.fq.gz # Built-in canonical encoding
./DBGPS2 -M examples/legacy-2base.map -a 1000 -b 1001 -o assembled reads.fq.gz
./DBGPS2 -M examples/3base-5bit.map -i 24 -l 120 -a 0 -b 1000 -o assembled reads.fq.gz
./DBGPS3 -M examples/3base-5bit.map -i 24 -l 120 -a 0 -b 1000 -o assembled reads.fq.gz
```

Use the codebook that was used to encode the input strands. In **DBGPS2**, the
encoded strand has this layout, with the same mapping applied to all sections:

```text
index | payload | CRC-16/XMODEM | zero padding to the next mapping block
```

The index is an unsigned integer serialized big-endian, including leading zero
bytes. CRC-16/XMODEM (polynomial `0x1021`, initial value zero, no reflection or
final XOR) covers those index bytes followed by the payload bytes. Its two bytes
are stored high byte first. Padding follows the CRC, is excluded from the CRC
calculation, and must be zero.

`-i` and `-l` remain lengths in **bases**. Both lengths must contain whole mapping
blocks and decode to whole bytes: `length % N == 0` and `(length / N * M) % 8 == 0`.
The CRC length is calculated automatically as `ceil(16/M) * N` bases. Thus, with
`3 bases -> 5 bits`, `-i 24` represents a 40-bit index, `-l 120` represents 200
payload bits, and the CRC occupies 12 bases with four trailing padding bits.

The index must fit within the k-mer (`1 <= k <= 31`), and `k - index_length`
must fit the existing 18-base forward primer. Index, payload, and encoded CRC
together must be shorter than 300 bases. `-a`/`-b` accept unsigned decimal values
within the index's bit capacity; use both options when the defaults do not fit.
An empty index (`-i 0 -a 0 -b 0`) is supported if the primer covers the k-mer.
For custom mappings the read counting limit adjusts to the full strand plus
the two existing 18-base primers.

**DBGPS3** searches the index plus the `-l` payload bases and returns their DNA
sequence. It validates codewords and requires a unique valid path, without CRC
screening or an automatically appended CRC region. Both index and payload must
contain complete codewords and whole decoded bytes. Its payload can be up to
300 bases, in addition to the index. If reads contain a following CRC, it is
outside this requested path length.

#### Mapping performance

The graph always stores `A=0, C=1, G=2, T=3`, preserving its rolling k-mer updates,
reverse complements, canonicalization, and coverage lookups. The configurable
codec is used only to construct index seeds and validate completed candidates.
Without `-M`, the index needs only a bit mask, and reconstructed candidate bytes
go directly to CRC validation. No mapping tables are allocated, loaded, or queried.

At startup, full-rate mappings with block sizes 1, 2, or 4 bases are compiled
into forward and reverse 256-entry byte tables. Custom candidate decoding therefore
uses one lookup per four bases. Other block sizes
use a dense DNA-word lookup table and a small integer bit buffer; no strings,
hash lookups, or codebook searches occur during candidate decoding. A forbidden
word rejects a candidate immediately. Temporary candidate buffers are on the
stack, eliminating two allocations per search that reaches CRC validation.

The general tables occupy `4 * 4^N + 2 * 2^M` bytes, plus metadata and the byte
tables; the maximum is 384 KiB. The 8-base limit bounds memory growth rather
than allocating tables exponential in an unrestricted block size. Runtime
depends mainly on graph branching and candidate count; the byte fast path does
not imply identical assembly times for different codebooks or code rates.

### DBGPS3

A next-generation de Bruijn graph assembler for DNA data storage with an interactive server mode.
The `-M` option uses the same codebooks as DBGPS2. Omitting it uses direct
`A=00,C=01,G=10,T=11` encoding; the original two-base permutation requires
`-M examples/legacy-2base.map`.

#### Command-Line Mode

```
DBGPS3 [options] <input file>
                  [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]

Options:
  -k INT     k-mer size [31]
  -i INT     length of index [16] bp
  -l INT     data encoding length [140] bp
  -t INT     number of threads [3]
  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [1]
  -C INT     Switch on k-mer coverage testing mode [0]
  -a INT     Initial index [101010102]
  -b INT     End index [101295684]
  -r FLOAT   Maximum coverage ratio [2.55]
  -R FLOAT   Upper bound for ratio range iteration (0 = disabled) [0.0]
  -s INT     Number of initial positions to skip ratio checks for [0]
  -o STR     Output prefix for log file and strands file
  -M FILE    Custom DNA/bit mapping (default: none; canonical A/C/G/T encoding)
```

#### Server Mode

DBGPS3 can also run as a TCP server, allowing remote clients to interactively query the de Bruijn graph using the JSON-RPC 2.0 protocol.

```bash
# Start server (default port 9300)
./DBGPS3 -S input.fq.gz

# Start on custom port with max 20 clients
./DBGPS3 -S 8080 -L 20 input.fq.gz

# Start with a custom codebook and compatible default query lengths
./DBGPS3 -S -M examples/3base-5bit.map -i 24 -l 120 input.fq.gz

# Connect with the interactive client
./dbgps3-client -h 127.0.0.1 -p 8080
```

For complete server mode documentation, including the JSON-RPC API reference, see [DBGPS3_SERVER.md](DBGPS3_SERVER.md).

**Server Options:**

| Flag | Description | Default |
|------|-------------|---------|
| `-S [PORT]` | Start in server mode | — |
| `-P INT` | Server TCP port | 9300 |
| `-L INT` | Max concurrent clients | 10 |

Query-time parameters (coverage cutoff, path length, index length, etc.) are provided per-request by the client, not at server startup. See [DBGPS3_SERVER.md](DBGPS3_SERVER.md) for details.

## License

This project is licensed under the terms in the LICENSE file.

## Author

Lifu Song (lifu.song@outlook.com)
