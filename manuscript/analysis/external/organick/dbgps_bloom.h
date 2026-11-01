/*
 * Blocked Bloom filter adapted from lh3/yak, bbf.c and yak-priv.h,
 * commit ef033d720e229593b32ad454af14ccea2120ec13.
 * https://github.com/lh3/yak
 *
 * The MIT License
 * Copyright (c) 2019-     Dana-Farber Cancer Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Local changes: checked allocation/sizing, automatic bits-per-key sizing,
 * full-width integer remixing, and a read-only membership operation. Unlike
 * The optional query filter preserves all counts. A separate two-pass
 * counting mode in dbgps_bloom_graph.h uses this primitive to drop singletons.
 */
#ifndef DBGPS_BLOOM_H
#define DBGPS_BLOOM_H

#include <stdint.h>
#include <stdlib.h>

#define DBGPS_BLOOM_BLOCK_SHIFT 9
#define DBGPS_BLOOM_BLOCK_MASK 511
#define DBGPS_BLOOM_HASHES 4
#define DBGPS_BLOOM_MAX_BITS_PER_KEY 64

typedef struct {
    int shift;
    size_t n_bytes;
    uint8_t *bits;
    void *allocation;
} dbgps_bloom_t;

/* The graph hash has only 2*k useful bits. Remix all 64 bits, including for
 * small k, before selecting a block and the two within-block hashes. */
static inline uint64_t dbgps_bloom_hash(uint64_t key)
{
    key = ~key + (key << 21);
    key ^= key >> 24;
    key = (key + (key << 3)) + (key << 8);
    key ^= key >> 14;
    key = (key + (key << 2)) + (key << 4);
    key ^= key >> 28;
    return key + (key << 31);
}

/* Strict CLI parsing: 0 disables; 1..64 specifies a minimum bits/key budget. */
static inline int dbgps_bloom_parse_bits(const char *text, int *bits)
{
    unsigned n = 0;
    if (!text || !*text || !bits) return -1;
    for (; *text; ++text) {
        if (*text < '0' || *text > '9') return -1;
        n = n * 10 + (unsigned)(*text - '0');
        if (n > DBGPS_BLOOM_MAX_BITS_PER_KEY) return -1;
    }
    *bits = (int)n;
    return 0;
}

static inline int dbgps_bloom_parse_shift(const char *text, int *shift)
{
    int n;
    if (dbgps_bloom_parse_bits(text, &n) < 0 || (n && (n < 23 || n > 40))) return -1;
    *shift = n;
    return 0;
}

static inline dbgps_bloom_t *dbgps_bloom_create(uint64_t n_keys, int bits_per_key)
{
    uint64_t needed;
    int shift = DBGPS_BLOOM_BLOCK_SHIFT;
    dbgps_bloom_t *b;
    /* Block selection plus the two 9-bit hashes must fit in a uint64_t. */
    if (bits_per_key < 1 || bits_per_key > DBGPS_BLOOM_MAX_BITS_PER_KEY ||
        n_keys > (UINT64_C(1) << 55) / (unsigned)bits_per_key) return NULL;
    needed = n_keys * (unsigned)bits_per_key;
    while ((UINT64_C(1) << shift) < needed) ++shift;
    if ((UINT64_C(1) << (shift - 3)) > SIZE_MAX - 63) return NULL;
    b = (dbgps_bloom_t*)calloc(1, sizeof(*b));
    if (!b) return NULL;
    b->shift = shift;
    b->n_bytes = (size_t)(UINT64_C(1) << (shift - 3));
    /* Preserve 64-byte alignment without platform-specific allocation APIs. */
    b->allocation = calloc(1, b->n_bytes + 63);
    if (!b->allocation) { free(b); return NULL; }
    b->bits = (uint8_t*)(((uintptr_t)b->allocation + 63) & ~(uintptr_t)63);
    return b;
}

static inline void dbgps_bloom_destroy(dbgps_bloom_t *b)
{
    if (!b) return;
    free(b->allocation);
    free(b);
}

/* A single builder writes the filter; after publication all queries are
 * read-only. Writers must finish before concurrent queries begin. */
static inline void dbgps_bloom_insert(dbgps_bloom_t *b, uint64_t key)
{
    uint64_t hash = dbgps_bloom_hash(key);
    int x = b->shift - DBGPS_BLOOM_BLOCK_SHIFT;
    uint64_t block = hash & ((UINT64_C(1) << x) - 1);
    unsigned z = (hash >> x) & DBGPS_BLOOM_BLOCK_MASK;
    unsigned h2 = (hash >> b->shift) & DBGPS_BLOOM_BLOCK_MASK;
    uint8_t *p = b->bits + (block << (DBGPS_BLOOM_BLOCK_SHIFT - 3));
    if ((h2 & 31) == 0) h2 = (h2 + 1) & DBGPS_BLOOM_BLOCK_MASK;
    for (int i = 0; i < DBGPS_BLOOM_HASHES; ++i) {
        p[z >> 3] |= (uint8_t)(1U << (z & 7));
        z = (z + h2) & DBGPS_BLOOM_BLOCK_MASK;
    }
}

/* NULL means disabled: every key proceeds to the exact hash table. */
static inline int dbgps_bloom_contains(const dbgps_bloom_t *b, uint64_t key)
{
    uint64_t hash, block;
    unsigned z, h2;
    int x;
    const uint8_t *p;
    if (!b) return 1;
    hash = dbgps_bloom_hash(key);
    x = b->shift - DBGPS_BLOOM_BLOCK_SHIFT;
    block = hash & ((UINT64_C(1) << x) - 1);
    z = (hash >> x) & DBGPS_BLOOM_BLOCK_MASK;
    h2 = (hash >> b->shift) & DBGPS_BLOOM_BLOCK_MASK;
    p = b->bits + (block << (DBGPS_BLOOM_BLOCK_SHIFT - 3));
    if ((h2 & 31) == 0) h2 = (h2 + 1) & DBGPS_BLOOM_BLOCK_MASK;
    for (int i = 0; i < DBGPS_BLOOM_HASHES; ++i) {
        if (!(p[z >> 3] & (1U << (z & 7)))) return 0;
        z = (z + h2) & DBGPS_BLOOM_BLOCK_MASK;
    }
    return 1;
}

#endif
