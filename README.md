# DBGPS2

DBGPS2 is a collection of tools for DNA data storage analysis and processing, with a focus on sequence assembly, cross-link detection, and coverage analysis.

## Overview

This project provides software tools for analyzing and reconstructing DNA sequences from data storage systems. The main components are:

- **DBGPS2**: A read-length constrained sequence assembler with dynamic coverage ratio support (formerly DBGPS-dy2-rl)
- **DBGPS-analyzer**: Tool for calculating key sequencing metrics (Sm, Kn, Kd values)
- **DBGPS-links**: Analyzer for detecting cross-links between DNA strands
- **DBGPS-seq-filter**: Filter to screen out entangled strands (also known as DBGPS-ft)
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
make DBGPS-analyzer
make DBGPS-links
make DBGPS-seq-filter
make DBGPS3
make dbgps3-server
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
```

### DBGPS-analyzer

Analyzer tool for calculating sequencing metrics like Sm, Kn, and Kd values.

#### Features
- Output prefix specification for log files
- Coverage ratio range iteration
- Parametrized read length and max coverage settings

#### Usage
```
DBGPS-analyzer [options] <Strand seq file> <NGS files>
                         [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]

Options:
  -k INT     k-mer size
  -t INT     number of threads
  -L INT     Maximal read length for k-mer counting
  -r FLOAT   Maximum coverage ratio (0 = no limitation)
  -R FLOAT   Upper bound for ratio range iteration
  -c INT     Minimum coverage cutoff
  -C INT     Maximum coverage cutoff
  -s INT     Number of initial ratio values to skip
  -o STR     Output prefix for additional files
```

### DBGPS-links

Analyzer for detecting cross-links between DNA strands in a dataset.

#### Usage
```
DBGPS-links [options] <in.fa>

Options:
  -k INT     k-mer size [31]
  -m INT     max link number [1]
```

### DBGPS-seq-filter

Filter tool to screen out entangled strands from DNA data. This tool identifies and removes strands with more than a specified number of cross-links, helping to improve sequence quality.

#### Features
- Primer-aware analysis that excludes primer regions from k-mer analysis
- Configurable k-mer size for entanglement detection
- Option to output either passed or filtered sequences

#### Usage
```
DBGPS-seq-filter [options] <in.fa>

Options:
  -k INT     k-mer size for entangle analysis [31]
  -m INT     Maximal strand cross-links [0]
  -p INT     Length of primers [18]
  -s         Output passed sequences instead of filtered ones
```

### DBGPS3

A next-generation de Bruijn graph assembler for DNA data storage with an interactive server mode.

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
```

#### Server Mode

DBGPS3 can also run as a TCP server, allowing remote clients to interactively query the de Bruijn graph using the JSON-RPC 2.0 protocol.

```bash
# Start server (default port 9300)
./DBGPS3 -S input.fq.gz

# Start on custom port with max 20 clients
./DBGPS3 -S 8080 -L 20 input.fq.gz

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
