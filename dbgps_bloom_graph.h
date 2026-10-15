/* Shared internals for the two assemblers' partitioned k-mer counters.
 * Include after kc_c4_t and DBGPS_COUNT_BITS/DBGPS_COUNT_MAX are defined. */
#ifndef DBGPS_BLOOM_GRAPH_H
#define DBGPS_BLOOM_GRAPH_H

#include <sys/stat.h>
#include "dbgps_bloom.h"

typedef struct {
    int p;
    kc_c4_t **h;
    dbgps_bloom_t *bloom; /* immutable query filter */
    dbgps_bloom_t **seen; /* partition-owned first-pass filters */
    int count_pass;       /* 0: exact; 1: candidates; 2: recount candidates */
} kc_c4x_t;

static void c4x_free_seen(kc_c4x_t *h)
{
    if (!h || !h->seen) return;
    for (int i = 0; i < 1 << h->p; ++i) dbgps_bloom_destroy(h->seen[i]);
    free(h->seen);
    h->seen = NULL;
}

static void c4x_destroy(kc_c4x_t *h)
{
    if (!h) return;
    c4x_free_seen(h);
    dbgps_bloom_destroy(h->bloom);
    if (h->h)
        for (int i = 0; i < 1 << h->p; ++i) kc_c4_destroy(h->h[i]);
    free(h->h);
    free(h);
}

static kc_c4x_t *c4x_init(int p)
{
    kc_c4x_t *h = (kc_c4x_t*)calloc(1, sizeof(*h));
    if (!h) return NULL;
    h->p = p;
    h->h = (kc_c4_t**)calloc((size_t)1 << p, sizeof(*h->h));
    if (!h->h) { c4x_destroy(h); return NULL; }
    for (int i = 0; i < 1 << p; ++i) {
        h->h[i] = kc_c4_init();
        if (!h->h[i]) { c4x_destroy(h); return NULL; }
    }
    return h;
}

/* Build only after exact counting; the query filter never removes a k-mer. */
static int c4x_set_bloom(kc_c4x_t *h, int bits_per_key)
{
    uint64_t n_keys = 0;
    dbgps_bloom_t *b = NULL;
    if (!h || bits_per_key < 0 || bits_per_key > DBGPS_BLOOM_MAX_BITS_PER_KEY) return -1;
    if (bits_per_key) {
        for (int i = 0; i < 1 << h->p; ++i) n_keys += kh_size(h->h[i]);
        b = dbgps_bloom_create(n_keys, bits_per_key);
        if (!b) return -1;
        for (int i = 0; i < 1 << h->p; ++i) {
            const kc_c4_t *table = h->h[i];
            for (khint_t j = 0; j < kh_end(table); ++j)
                if (kh_exist(table, j))
                    dbgps_bloom_insert(b, ((kh_key(table, j) >> DBGPS_COUNT_BITS) << h->p) | (uint64_t)i);
        }
    }
    dbgps_bloom_destroy(h->bloom);
    h->bloom = b;
    return 0;
}

/* Called by the existing partition workers: each filter has one writer. */
static inline void c4x_count_key(kc_c4x_t *h, int partition, uint64_t hash)
{
    kc_c4_t *table = h->h[partition];
    uint64_t key = hash >> h->p << DBGPS_COUNT_BITS;
    khint_t pos;
    int absent;
    if (h->count_pass == 1) {
        int seen = dbgps_bloom_contains(h->seen[partition], hash);
        dbgps_bloom_insert(h->seen[partition], hash);
        if (seen) kc_c4_put(table, key, &absent); /* count remains zero */
        return;
    }
    if (h->count_pass == 2) {
        pos = kc_c4_get(table, key);
        if (pos == kh_end(table)) return;
    } else pos = kc_c4_put(table, key, &absent);
    if ((kh_key(table, pos) & DBGPS_COUNT_MAX) < DBGPS_COUNT_MAX) ++kh_key(table, pos);
}

static kc_c4x_t *count_file2(const char *fn, kc_c4x_t *h, int k, int p,
                            int block_size, int n_thread, int read_len);

static void c4x_prune_singletons(kc_c4x_t *h)
{
    for (int i = 0; i < 1 << h->p; ++i) {
        kc_c4_t *table = h->h[i];
        for (khint_t j = 0; j < kh_end(table);) {
            if (kh_exist(table, j) && (kh_key(table, j) & DBGPS_COUNT_MAX) < 2)
                kc_c4_del(table, j); /* deletion may move another key here */
            else ++j;
        }
    }
}

/* -F is log2 of the TOTAL bitmap bits across partitions. All input files
 * participate in both passes, including repeats spanning file boundaries. */
static kc_c4x_t *count_files_filtered(int n_files, char **files, int k, int p,
                                     int block_size, int n_thread, int read_len, int shift)
{
    kc_c4x_t *h;
    if (shift < p + DBGPS_BLOOM_BLOCK_SHIFT || shift > 40) return NULL;
    /* Two-pass counting requires replayable inputs. Reject FIFOs/devices
     * before opening them, so a second pass cannot block or silently empty. */
    for (int i = 0; i < n_files; ++i) {
        struct stat st;
        if (!files[i] || stat(files[i], &st) < 0 || !S_ISREG(st.st_mode)) return NULL;
    }
    h = c4x_init(p);
    if (!h) return NULL;
    h->seen = (dbgps_bloom_t**)calloc((size_t)1 << p, sizeof(*h->seen));
    if (!h->seen) { c4x_destroy(h); return NULL; }
    for (int i = 0; i < 1 << p; ++i) {
        h->seen[i] = dbgps_bloom_create(UINT64_C(1) << (shift - p), 1);
        if (!h->seen[i]) { c4x_destroy(h); return NULL; }
    }
    for (int pass = 1; pass <= 2; ++pass) {
        h->count_pass = pass;
        fprintf(stderr, "Bloom counting pass %d/2 (%s) ......\n", pass,
                pass == 1 ? "select candidates" : "exact recount");
        for (int i = 0; i < n_files; ++i)
            if (!count_file2(files[i], h, k, p, block_size, n_thread, read_len)) {
                c4x_destroy(h);
                return NULL;
            }
        if (pass == 1) c4x_free_seen(h);
    }
    h->count_pass = 0;
    /* Bloom false positives can admit singletons as candidates. Removing
     * exact counts < 2 makes the final graph independent of those errors. */
    c4x_prune_singletons(h);
    return h;
}

#endif
