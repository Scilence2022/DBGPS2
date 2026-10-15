# Redundancy and error correction

DBGPS2 and DBGPS3 share the same byte-level protection module and the same
`-E CODE` option. Both default to `crc16`. For DBGPS3 inputs without a checksum,
explicitly select `-E none` to retain unprotected path search. No new runtime
library is required.

| `-E` value | Trailer | Behavior |
|------------|---------|----------|
| `none` | None | Codebook validation and path uniqueness only |
| `crc16` | 2 bytes | CRC-16/XMODEM detection, compatible with existing DBGPS2 data |
| `crc32` | 4 bytes | CRC-32/ISO-HDLC (IEEE) detection |
| `crc32c` | 4 bytes | CRC-32/ISCSI (Castagnoli) detection |
| `rs:NB` (or `rs`) | N bytes (`rs` means `rs:16B`) | Shortened systematic Reed–Solomon; corrects up to `floor(N/2)` erroneous bytes |
| `hamming` | `ceil(5D/8)` bytes for D data bytes | Extended Hamming(13,8) SECDED; corrects one bit and detects two bits per protected data byte plus its five parity bits |

Examples (the input must have been encoded with the selected protection):

```bash
./DBGPS2 -E rs:8B -a 1000 -b 1001 -o assembled reads.fa
./DBGPS2 -E crc32c -a 1000 -b 1001 -o assembled reads.fa
./DBGPS3 -E rs:8B -a 1000 -b 1001 -o assembled reads.fa
./DBGPS3 -E none -a 1000 -b 1001 -o assembled reads.fa
./DBGPS3 -S -E hamming reads.fa
./DBGPS3 -S -E rs:4B -M examples/3base-5bit.map -i 24 -l 120 reads.fa
```

Use uppercase `B` to specify bytes: `rs:8B` means eight parity bytes (64 bits).
The legacy `rs:8` spelling is accepted with the same meaning. Lowercase `b`,
including `rs:8b`, is rejected. Bare `rs` selects `rs:16B`. Help text, CLI logs,
and the server's `redundancy.spec` field use the explicit byte-unit form.
CRC names remain standard algorithm names: the 16 in `crc16` denotes bits.

## Wire format and lengths

All protection is systematic and covers **index bytes followed by payload
bytes**, before DNA mapping. The index is an unsigned big-endian integer,
including leading zero bytes. With any codebook, the layout is:

```text
index | payload | redundancy trailer | zero bits to fill the last mapping block
```

`-i` and `-l` (JSON-RPC `indexLen` and `pathLen`) describe index and payload
lengths in bases; neither includes redundancy. Both must consist of complete
mapping blocks and decode to whole bytes. For a mapping of B bases to M bits,
P trailer bytes require `ceil(8P/M) * B` extra bases. The final mapping padding
must be zero and is excluded from protection. Hamming also has its own zero
padding inside the trailer, which is validated separately.

DBGPS2 requires `index + payload + redundancy` to be shorter than 300 bases.
DBGPS3 permits `payload + redundancy` up to 300 bases, plus the index. Both
validate lengths and codeword capacities before counting reads. The read
counting limit grows to include the complete strand and both 18-base primers
when necessary, including canonical mappings.

The output contains the index, payload, trailer and mapping padding as DNA.
For `none`, it contains only the index and payload. Error-correcting codes
re-encode corrected bytes using the selected codebook, including repaired
parity. Multiple graph paths that produce identical corrected strands count
as one result; different valid corrected strands are ambiguous and rejected.
A correction that changes the seeded index is rejected.

ECC operates on completed graph candidates. It cannot recover a missing seed
or a path absent from the graph, and it does not itself correct insertions or
deletions. Forbidden codebook words and nonzero mapping padding fail before
ECC. Beyond their stated correction/detection bounds, RS and Hamming can
miscorrect or accept another valid codeword; they do not provide an additional
independent CRC check. These are inner codes per strand, not outer codes
across strands or erasure decoders.

## Exact code parameters

### Checksums

Checksum bytes are stored most significant byte first, independent of the
reflection used during CRC computation.

| Code | Polynomial | Initial value | Reflected input/output | Final XOR | Check for `123456789` |
|------|------------|---------------|------------------------|-----------|-----------------------|
| `crc16` | `0x1021` | `0x0000` | No | `0x0000` | `0x31c3` |
| `crc32` | `0x04c11db7` | `0xffffffff` | Yes | `0xffffffff` | `0xcbf43926` |
| `crc32c` | `0x1edc6f41` | `0xffffffff` | Yes | `0xffffffff` | `0xe3069283` |

### Reed–Solomon

- GF(256), primitive polynomial `0x11d`, primitive element `alpha = 2`.
- Generator `g(x) = product(x + alpha^i, i=0..N-1)`; first consecutive root 0.
- Bytes are polynomial coefficients, highest degree first. Parity follows data.
- `1 <= N <= 254`, at least one data byte, and `data_bytes + N <= 255`.
- A shortened codeword is equivalent to omitting leading zero data symbols.
  No implicit chunking or interleaving is performed.
- Decoding uses Berlekamp–Massey, Chien search over the shortened word,
  Forney error magnitudes, and a final syndrome check.

The format matches
[`reedsolo.RSCodec(N, nsize=255, fcr=0, prim=0x11d, generator=2)`](https://github.com/tomerfiliba-org/reedsolomon).
Other RS encoders must use these exact parameters. For example, `hello world`
with N=10 has trailer `ed2554c4fdfd89f3a8aa`.

### Hamming SECDED

Each data byte is an independent extended Hamming(13,8) word. Place its bits
MSB first at positions `3,5,6,7,9,10,11,12`. Even parity occupies positions
`1,2,4,8`; position 13 is overall even parity. Keep the original data bytes
contiguous and append their parity groups in data-byte order. Each five-bit
group is `p8,p4,p2,p1,p13`, packed MSB first into the trailer; pad its final
byte with zero bits. This systematic trailer packing is part of DBGPS's wire
format and must match the writer. A single error in any of the thirteen bits
is corrected; any double error within that word is rejected.

## DBGPS3 server and C API

`-E` selects immutable protection shared by all queries. `graph/info` exposes
`redundancy.code`, `spec`, `parity_bytes`, `parity_bases`, `padding_bits`, and
`path_bases` for the startup layout. `path_bases` includes payload and trailer;
the existing top-level `path_len` remains the payload length. `code` remains
the algorithm name (`rs`), while `spec` includes its parity size (`rs:8B`).

`path/find` uses the selected protection with each request's index and payload
lengths. `path/findKmer` still returns only the extension after its starting
k-mer. With protection enabled, its last `indexLen` bases are interpreted as
the index (default: startup `-i`), included in validation, then omitted from
the response. Set `indexLen: 0` for a format whose protection covers only the
extension. With `none`, the original arbitrary-k-mer behavior is unchanged.
`graph/load` preserves both mapping and redundancy configuration.

For C callers, configure `g->index_len` and `g->path_len`, then call
`dbgps3_set_mapping()` and `dbgps3_set_redundancy()` before serving queries.
Both setters preserve the previous configuration on failure. A zero-initialized
`dbgps_redundancy_t` selects `none`; graphs created by `dbgps3_load()` or
`dbgps3_load_filtered()` explicitly initialize protection to `crc16`. The byte API in `dbgps_redundancy.h` can
also encode input data and validate/correct complete codewords independently
of the graph. Decode failures leave its input unchanged. DNA encoding and
decoding use `dbgps_codec.h`, and `dbgps_strand.h` composes both layers.

## Adding another code

Implement three byte-level callbacks in `dbgps_redundancy.c`: calculate trailer
size/capacity, append parity, and validate/correct in place. Add a descriptor
to its registry and extend the option help and parameter parser if needed.
Callbacks must be reentrant, avoid mutable shared state, validate before
modifying input on failure, and produce the documented systematic format.
No graph traversal or DNA mapping changes are needed. Add independent format
vectors, boundary and failure tests, and exercise both assemblers and the
server. Codes requiring another wire layout or cross-strand state need an
explicit extension to this interface.

`make test` includes checksum vectors, an independent RS fixture, all 254 RS
parity lengths, every Hamming byte value with every single/double bit error,
independent Python encoders, index protection, corrected-result uniqueness,
all 72 DNA mapping widths, CLI correction and server reloads/concurrent queries.
Use `make test asan=1` for AddressSanitizer.
