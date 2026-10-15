# DBGPS3 Server Mode

DBGPS3 supports a **server mode** that exposes the de Bruijn graph engine over a TCP interface using the **JSON-RPC 2.0** protocol. This allows remote clients to interactively query k-mer coverage, find paths, load additional data, and retrieve graph information without restarting the process.

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                       DBGPS3 Server                      │
│                                                          │
│  ┌─────────────────────────────────────────────────────┐ │
│  │              dbgps3_core (Graph Engine)              │ │
│  │                                                      │ │
│  │  • dbgps3_load()          – Build de Bruijn graph   │ │
│  │  • dbgps3_kmer_cov()      – Query k-mer coverage    │ │
│  │  • dbgps3_find_path()     – Decode DNA path         │ │
│  │  • dbgps3_load_more()     – Hot-load more data      │ │
│  │  • dbgps3_destroy()       – Free all memory         │ │
│  └─────────────────────────────────────────────────────┘ │
│                          ▲                               │
│                          │ thread-safe                   │
│                          ▼                               │
│  ┌─────────────────────────────────────────────────────┐ │
│  │           dbgps3_server (JSON-RPC 2.0)               │ │
│  │                                                      │ │
│  │  • Multi-client support (configurable)              │ │
│  │  • One thread per client                            │ │
│  │  • Graceful shutdown via SIGINT / RPC               │ │
│  └─────────────────────────────────────────────────────┘ │
│                          ▲                               │
│                          │ TCP (port 9300 default)       │
└──────────────────────────┼───────────────────────────────┘
                           │
              ┌────────────┴────────────┐
              │     dbgps3-client       │
              │   (interactive CLI)     │
              └─────────────────────────┘
```

### Key Design Decisions

| Aspect | Design |
|--------|--------|
| **No global state** | All graph data is in `dbgps3_graph_t *`, enabling concurrent client access |
| **Thread safety** | `find_paths` allocates per-call buffers; the hash table is read-only after build |
| **Output abstraction** | Results are written into a caller-provided buffer, not directly to `stdout`/`stderr` |
| **JSON-RPC 2.0** | Standard protocol with `id`, `method`, `params`, `result`, `error` fields |
| **Hot reload** | `graph/load` command allows adding data without server restart |

## Starting the Server

### Build

```bash
make DBGPS3
```

### Launch

```bash
# Basic server (default port 9300, CRC16 protection)
./DBGPS3 -S input.fq.gz

# Input without a checksum
./DBGPS3 -S -E none input.fq.gz

# Specify port
./DBGPS3 -S 8080 input.fq.gz
# or
./DBGPS3 -S -P 8080 input.fq.gz

# Max concurrent clients (default: 10)
./DBGPS3 -S -P 8080 -L 20 input.fq.gz

# With graph build options
./DBGPS3 -S -k 31 -t 4 input.fq.gz

# Custom mapping shared by every query (same format as DBGPS2)
./DBGPS3 -S -M examples/3base-5bit.map -i 24 -l 120 input.fq.gz

# Enable Reed–Solomon correction (input must use matching RS parameters)
./DBGPS3 -S -E rs:8B input.fq.gz
```

### Server-Specific Options

| Flag | Description | Default |
|------|-------------|---------|
| `-S [PORT]` | Start server mode. Optional port argument. | 9300 |
| `-P INT` | Server TCP port. Overrides `-S` port. | 9300 |
| `-L INT` | Max concurrent clients. | 10 |
| `-B INT` | Query Bloom bitmap bits per retained k-mer (1..64; 0 disables). | 0 |
| `-F INT` | Two-pass singleton filtering; log2 total bitmap bits (23..40; 0 disables). | 0 |
| `-M FILE` | Immutable DNA/bit codebook, shared by all queries. | No custom mapping |
| `-E CODE` | `none`, `crc16`, `crc32`, `crc32c`, `rs[:NB]`, `hamming`; shared by all queries. | `crc16` |

### Graph Build Options

These options affect the de Bruijn graph construction and are set at server startup:

| Flag | Description | Default |
|------|-------------|---------|
| `-k INT` | K-mer size. | 31 |
| `-t INT` | Number of threads for graph building. | 3 |
| `-p INT` | Prefix length in bits (must be 14). | 14 |

Without `-M`, index seeds use the canonical k-mer encoding directly:
`A=00,C=01,G=10,T=11`. Use `-M examples/legacy-2base.map` for data encoded with
the original special two-base permutation. The codebook is parsed once before
graph construction and remains read-only during queries. There is no per-query
codebook filename parameter.

All DBGPS2 codebook widths are supported: 1–8 bases to 1–16 bits, with at most
two bits per base. Index and payload lengths must contain complete mapping
blocks and decode to whole bytes. For a 3-base/5-bit code, `-i 24` is a 40-bit
index and `-l 120` is a 200-bit payload. By default (`-E crc16`), DBGPS3
extends the search to include the CRC16 trailer and checks index plus payload
bytes. Explicit `-E none` searches only the index and payload without checksum
screening. Other `-E` modes use their corresponding trailers and validation/correction. `pathLen` always excludes redundancy. The payload plus trailer must
fit within 300 bases. Different candidates corrected to the same DNA count as
one result; distinct corrected strands are rejected as ambiguous. See
[REDUNDANCY.md](REDUNDANCY.md) for exact wire formats and correction limits.

### Per-Query Parameters

Clients may provide these parameters in each request. Startup `-i` and `-l`
set the default index and payload lengths for requests that omit them:

| Parameter | API Key | Description | Default |
|-----------|---------|-------------|---------|
| Coverage cutoff | `cutCov` | Minimum k-mer coverage to follow. | 1 |
| Max coverage ratio | `maxCovRatio` | Max ratio for branching decisions. | 2.55 |
| Skip ratios | `skipRatios` | Number of initial positions to skip ratio checks. | 0 |
| Payload length | `pathLen` | Payload bases; the selected redundancy is added automatically. | 140 |
| Index length | `indexLen` | Index bases; also used by protected `path/findKmer`. | 16 |

## Client Usage

### Build

```bash
make dbgps3-client
```

### Interactive Mode

```bash
./dbgps3-client
# or specify host/port
./dbgps3-client -h 192.168.1.100 -p 8080
```

```
dbgps3> graph/info
[Response...]

dbgps3> kmer/cov ACGTACGT...
[Response...]

dbgps3> find 101010102
[Response...]

dbgps3> exit
```

### Shortcut Command Mode

```bash
# One-shot shortcut command
./dbgps3-client -c 'graph/info'
./dbgps3-client -c 'kmer/cov ACGTACGTACGTACGTACGTACGT'
./dbgps3-client -c 'find 101010102'
./dbgps3-client -c 'findKmer 101010102 5 2.5 0'

# Pipe raw JSON
echo '{"jsonrpc":"2.0","method":"graph/info","id":1}' | ./dbgps3-client
```

### Client Options

| Flag | Description | Default |
|------|-------------|---------|
| `-h HOST` | Server host. | 127.0.0.1 |
| `-p PORT` | Server port. | 9300 |
| `-c 'CMD'` | Shortcut command (one-shot mode). | — |

## JSON-RPC 2.0 API Reference

All requests follow the JSON-RPC 2.0 specification. Requests with `"id": null` are **notifications** (no response expected).

### 1. `kmer/cov` — Query K-mer Coverage

Get the coverage count for a k-mer in the graph.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "kmer/cov",
  "params": {
    "kmer": "ACGTACGTACGTACGTACGTACGT"
  },
  "id": 1
}
```

**Success Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "kmer": "ACGTACGTACGTACGTACGTACGT",
    "key": 1234567890123456789,
    "coverage": 42
  },
  "id": 1
}
```

**K-mer Not Found:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "kmer": "ACGTACGTACGTACGTACGTACGT",
    "coverage": 0
  },
  "id": 1
}
```

### 2. `path/find` — Find Path by Index

Decode a DNA strand by its decimal index number. The server encodes the index
into a starting k-mer using `indexLen` and the startup codebook. `index` accepts
an unsigned integer up to `2^53-1` or an unsigned decimal string for exact larger
values. The index must fit the selected layout. The client `find <index>`
shortcut sends a decimal string, preserving all supported index bits.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "path/find",
  "params": {
    "index": 101010102,
    "cutCov": 1,
    "maxCovRatio": 2.55,
    "skipRatios": 0,
    "indexLen": 16,
    "pathLen": 140
  },
  "id": 2
}
```

All parameters except `index` are optional and default to the graph's built-in values.

**Success Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "index": 101010102,
    "found": true,
    "sequence": "ACGTACGT..."
  },
  "id": 2
}
```

**Path Not Found:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "index": 101010102,
    "found": false
  },
  "id": 2
}
```

String indices are returned as decimal strings. Incompatible block lengths,
out-of-range indices, fractional numbers, and invalid query parameters return
JSON-RPC error `-32602` before constructing a primer seed.

### 3. `path/findKmer` — Find Path by Raw K-mer

Find a path starting from an arbitrary k-mer sequence. Since the client provides
the starting k-mer directly, `indexLen` is not needed with `-E none`. The literal
k-mer always uses canonical A/C/G/T encoding; the returned extension must start at a codeword
boundary and satisfy the startup codebook. Only the extension is returned.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "path/findKmer",
  "params": {
    "kmer": "ACGTACGTACGTACGTACGTACGTACGTACG",
    "cutCov": 1,
    "maxCovRatio": 2.55,
    "skipRatios": 0,
    "pathLen": 140
  },
  "id": 3
}
```

All parameters except `kmer` are optional. The k-mer must contain exactly `k` bases.

**Success Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "kmer": "ACGTACGTACGTACGTACGTACGTACGTACG",
    "found": true,
    "sequence": "ACGTACGT..."
  },
  "id": 3
}
```

With protection enabled, the last `indexLen` bases of the k-mer are included
as index bytes in validation/correction (default: startup `-i`). The returned
sequence still excludes that index and contains the payload plus redundancy.
Use `indexLen: 0` when the encoded protection covers only the extension.

### 4. `graph/info` — Get Graph Information

The response also includes `redundancy`: `code`, `spec`, `parity_bytes`, `parity_bases`,
`padding_bits`, and `path_bases` for the default query layout. The existing
`path_len` is the payload length; `redundancy.path_bases` includes its trailer.
For RS, `code` is `rs` and `spec` is the canonical byte-unit form, e.g.
`rs:8B`. Legacy `rs:8` startup arguments remain accepted; lowercase `rs:8b`
is rejected. Mapping and protection remain unchanged across `graph/load`.

Retrieve graph metadata.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "graph/info",
  "params": {},
  "id": 4
}
```

**Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "k": 31,
    "kmer_count": 12345678,
    "index_len": 16,
    "index_byte_num": 4,
    "path_len": 140,
    "bloom": {
      "enabled": false,
      "bits_per_key": 0,
      "bitmap_bytes": 0,
      "count_shift": 0,
      "singletons_removed": false
    },
    "mapping": {
      "custom": false,
      "bases": 1,
      "bits": 2,
      "mode": "canonical"
    },
    "redundancy": {
      "code": "crc16",
      "spec": "crc16",
      "parity_bytes": 2,
      "parity_bases": 8,
      "padding_bits": 0,
      "path_bases": 148
    }
  },
  "id": 4
}
```

`index_len` and `path_len` are the default query lengths. `index_byte_num` counts
decoded data bytes; for a 24-base/40-bit index it is 5. Mapping mode is
`canonical`, `byte_lookup`, or `block_lookup`. Clients may override lengths
per query, subject to the same codeword and byte alignment checks.

`bloom.enabled` reports the query filter controlled by `-B`; positive probes
still use exact coverage lookup. `count_shift` independently reports the `-F`
two-pass mode, which removes singletons. Both are disabled by default. See
[BLOOM_FILTER.md](BLOOM_FILTER.md) for configuration, implementation selection,
and performance measurements.

### 5. `graph/load` — Load Additional Data

Hot-load more FASTQ/FASTA files into the existing graph.

Loads hold exclusive graph access while concurrent queries share read access.
With `-B`, the query filter is rebuilt after successful loading. With `-F`, the
server rebuilds from all historical inputs plus the new files so a prior
singleton can become a correctly counted repeated k-mer. Keep those regular
input files available and unchanged. A failed filtered rebuild preserves the
old graph; reload temporarily holds both the old and replacement graphs.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "graph/load",
  "params": {
    "files": ["/path/to/new_data.fq.gz", "/path/to/more_data.fa"]
  },
  "id": 5
}
```

**Success Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "status": "ok",
    "message": "Successfully loaded 2 files"
  },
  "id": 5
}
```

**Error Response:**
```json
{
  "jsonrpc": "2.0",
  "error": {
    "code": -32603,
    "message": "Internal error",
    "data": "Failed to open file: /path/to/missing.fq"
  },
  "id": 5
}
```

### 6. `server/shutdown` — Shutdown Server

Gracefully stop the server.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "server/shutdown",
  "params": {},
  "id": 6
}
```

**Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "status": "ok",
    "message": "Server shutting down"
  },
  "id": 6
}
```

### Error Responses

All errors follow the JSON-RPC 2.0 error format:

```json
{
  "jsonrpc": "2.0",
  "error": {
    "code": -32600,
    "message": "Invalid Request",
    "data": "missing method field"
  },
  "id": null
}
```

| Code | Message | Meaning |
|------|---------|---------|
| -32600 | Invalid Request | Malformed JSON or missing fields |
| -32601 | Method not found | Unknown method name |
| -32602 | Invalid params | Invalid parameter structure |
| -32603 | Internal error | Server-side processing error |

## Building a Custom Client

Any language that supports TCP sockets can be used as a client. Here's a Python example:

```python
import socket
import json

def dbgps3_request(method, params=None, host="127.0.0.1", port=9300):
    request = {
        "jsonrpc": "2.0",
        "method": method,
        "params": params or {},
        "id": 1
    }

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.connect((host, port))
        s.sendall((json.dumps(request) + "\n").encode())
        response = b""
        while True:
            chunk = s.recv(4096)
            if not chunk:
                break
            response += chunk
        return json.loads(response.decode())

# Example usage
result = dbgps3_request("graph/info")
print(f"K-mer count: {result['result']['kmer_count']}")

result = dbgps3_request("kmer/cov", {"kmer": "ACGT" * 8})
print(f"Coverage: {result['result']['coverage']}")
```

## Programmatic C API

For embedding the DBGPS3 engine into other C applications:

```c
#include "dbgps3_core.h"
#include "dbgps3_server.h"
#include <stdio.h>

/* Build graph (CRC16 protection by default; configure "none" explicitly
   with dbgps3_set_redundancy() for input without a checksum). */
dbgps3_graph_t *g = dbgps3_load(31, 14, 3, 10000000, 200, 1, &input_file);
if (!g) return 1;

/* Optional custom mapping: set compatible defaults before configuring it.
   The graph owns the codebook; configure it before starting any query threads. */
g->index_len = 24;
g->path_len = 120;
char error[512];
if (dbgps3_set_mapping(g, "examples/3base-5bit.map", error, sizeof(error)) < 0) {
    fprintf(stderr, "%s\n", error);
    dbgps3_destroy(g);
    return 1;
}

uint64_t encoded_index;
dbgps3_encode_index(g, 1000, g->index_len, &encoded_index);
uint64_t start_kmer = dbgps3_primer_uint64("CCTGCAGAGTAGCATGTC", 18, g->index_len, g->k)
                      | encoded_index;

/* Query k-mer coverage */
uint64_t key = dbgps3_actgkmer_hashkey("ACGTACGTACGTACGTACGTACGTACGTACG", 31);
unsigned short cov = dbgps3_kmer_cov(g, key);

/* Find a path (path_len and index_byte_num are explicit per-query params) */
char buf[DBGPS3_MAX_PATH_LEN + 128];
int rc = dbgps3_find_path(g, start_kmer, 1, 2.55, 0,
                          g->path_len, g->index_byte_num,
                          buf, sizeof(buf));

/* Start server */
dbgps3_server_t srv = { .graph = g, .port = 9300, .max_clients = 10 };
dbgps3_server_start(&srv);

/* Cleanup */
dbgps3_destroy(g);
```

## Performance Considerations

- **Concurrent queries**: Queries share read access; `graph/load` holds a write lock while updating the graph and its Bloom filter.
- **Memory**: The graph is shared across all clients. Monitor memory usage with large datasets.
- **Network**: Responses for `path/find` can be large (up to ~300 bytes for the DNA sequence). The protocol uses `\n` as the message delimiter.
- **Client limit**: The `-L` flag controls maximum concurrent connections. Excess connections are rejected.
