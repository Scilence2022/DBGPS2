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
make dbgps3-server
```

### Launch

```bash
# Basic server (default port 9300)
./DBGPS3 -S input.fq.gz

# Specify port
./DBGPS3 -S 8080 input.fq.gz
# or
./DBGPS3 -S -P 8080 input.fq.gz

# Max concurrent clients (default: 10)
./DBGPS3 -S -P 8080 -L 20 input.fq.gz

# With graph build options
./DBGPS3 -S -k 31 -t 4 input.fq.gz
```

### Server-Specific Options

| Flag | Description | Default |
|------|-------------|---------|
| `-S [PORT]` | Start server mode. Optional port argument. | 9300 |
| `-P INT` | Server TCP port. Overrides `-S` port. | 9300 |
| `-L INT` | Max concurrent clients. | 10 |

### Graph Build Options

These options affect the de Bruijn graph construction and are set at server startup:

| Flag | Description | Default |
|------|-------------|---------|
| `-k INT` | K-mer size. | 31 |
| `-t INT` | Number of threads for graph building. | 3 |
| `-p INT` | Prefix length in bits (min: 14). | 14 |

### Per-Query Parameters

The following parameters are provided by the client in each request, not at server startup:

| Parameter | API Key | Description | Default |
|-----------|---------|-------------|---------|
| Coverage cutoff | `cutCov` | Minimum k-mer coverage to follow. | 1 |
| Max coverage ratio | `maxCovRatio` | Max ratio for branching decisions. | 2.55 |
| Skip ratios | `skipRatios` | Number of initial positions to skip ratio checks. | 0 |
| Path length | `pathLen` | How far to extend the path search (in bases). | 140 |
| Index length | `indexLen` | Index encoding length (for `path/find` only). | 16 |

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

Decode a DNA strand by its decimal index number. The server encodes the index into a starting k-mer using `indexLen`.

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
    "path": "CCTGCAGAGTAGCATGTCACGTACGT...",
    "c_index": 42
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
    "path": "",
    "c_index": -1
  },
  "id": 2
}
```

### 3. `path/findKmer` — Find Path by Raw K-mer

Find a path starting from an arbitrary k-mer sequence. Since the client provides the starting k-mer directly, `indexLen` is not needed.

**Request:**
```json
{
  "jsonrpc": "2.0",
  "method": "path/findKmer",
  "params": {
    "startKmer": "ACGTACGTACGTACGTACGTACGT",
    "cutCov": 1,
    "maxCovRatio": 2.55,
    "skipRatios": 0,
    "pathLen": 140
  },
  "id": 3
}
```

All parameters except `startKmer` are optional.

**Success Response:**
```json
{
  "jsonrpc": "2.0",
  "result": {
    "startKmer": "ACGTACGTACGTACGTACGTACGT",
    "path": "ACGTACGTACGTACGTACGTACGT...",
    "c_index": 42
  },
  "id": 3
}
```

### 4. `graph/info` — Get Graph Information

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
    "kmerSize": 31,
    "prefixBits": 14,
    "kmerCount": 12345678,
    "indexLen": 16,
    "indexByteNum": 4,
    "pathLen": 140,
    "threadCount": 3,
    "blockSize": 10000000,
    "readLen": 200
  },
  "id": 4
}
```

> **Note:** `indexLen` and `pathLen` in `graph/info` are the graph's default values. Clients can override them per-query in `path/find` and `path/findKmer`.

### 5. `graph/load` — Load Additional Data

Hot-load more FASTQ/FASTA files into the existing graph.

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
print(f"K-mer count: {result['result']['kmerCount']}")

result = dbgps3_request("kmer/cov", {"kmer": "ACGT" * 8})
print(f"Coverage: {result['result']['coverage']}")
```

## Programmatic C API

For embedding the DBGPS3 engine into other C applications:

```c
#include "dbgps3_core.h"
#include "dbgps3_server.h"

/* Build graph */
dbgps3_graph_t *g = dbgps3_load(31, 14, 3, 10000000, 200, 1, &input_file);

/* Query k-mer coverage */
uint64_t key = dbgps3_actgkmer_hashkey("ACGTACGTACGTACGTACGTACGT", 31);
unsigned short cov = dbgps3_kmer_cov(g, key);

/* Find a path (path_len and index_byte_num are explicit per-query params) */
char buf[DBGPS3_MAX_PATH_LEN + 128];
int rc = dbgps3_find_path(g, start_kmer, 1, 2.55, 0,
                          140, 4,  /* path_len, index_byte_num */
                          buf, sizeof(buf));

/* Start server */
dbgps3_server_t srv = { .graph = g, .port = 9300, .max_clients = 10 };
dbgps3_server_start(&srv);

/* Cleanup */
dbgps3_destroy(g);
```

## Performance Considerations

- **Concurrent queries**: The hash table is read-only after initial build, so multiple clients can query simultaneously without locking.
- **Memory**: The graph is shared across all clients. Monitor memory usage with large datasets.
- **Network**: Responses for `path/find` can be large (up to ~300 bytes for the DNA sequence). The protocol uses `\n` as the message delimiter.
- **Client limit**: The `-L` flag controls maximum concurrent connections. Excess connections are rejected.
