/* Standalone synthetic benchmark; see BLOOM_FILTER.md for commands and limits. */
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "dbgps_core.h"

static uint64_t key_for(uint64_t i)
{
    /* Distinct packed words with scattered prefixes and reverse complements. */
    i *= UINT64_C(0x9e3779b97f4a7c15);
    return dbgps_min_hash_key(i & ((UINT64_C(1) << 62) - 1), 31);
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    if (argc == 3 && !strcmp(argv[1], "generate")) {
        FILE *out = fopen(argv[2], "w");
        assert(out);
        char dna[32];
        dna[31] = 0;
        for (unsigned i = 0; i < 1000000; ++i) {
            dbgps_uint64_acgt(key_for(i), dna, 31);
            int copies = i < 100000 ? 3 : 1;
            for (int j = 0; j < copies; ++j) fprintf(out, ">%u_%d\n%s\n", i, j, dna);
        }
        fclose(out);
        return 0;
    }
    if (argc != 4) {
        fprintf(stderr, "Usage: bench_bloom generate reads.fa | reads.fa QUERY_BITS COUNT_SHIFT\n");
        return 1;
    }
    int bits = atoi(argv[2]), shift = atoi(argv[3]);
    double start = now();
    dbgps_graph_t *g = dbgps_load_filtered(31, DBGPS_KC_BITS, 3, 100000, 31, 1, argv + 1, shift);
    assert(g && dbgps_set_bloom(g, bits) == 0);
    double build = now() - start;
    uint64_t sum = 0;
    const unsigned queries = 10000000;
    start = now();
    for (unsigned i = 0; i < queries; ++i) {
        /* 5% repeated present keys, 95% absent. Generate the same queries in
         * all modes; the checksum must agree despite singleton filtering. */
        uint64_t n = i % 20 == 0 ? i % 100000 : 1000000 + (uint64_t)i;
        sum += dbgps_kmer_cov(g, key_for(n));
    }
    printf("bits=%d shift=%d kmers=%ld bitmap=%zu build_s=%.3f query_ns=%.2f checksum=%"PRIu64"\n",
            bits, shift, dbgps_get_kmer_count(g), dbgps_bloom_bytes(g), build,
            (now() - start) * 1e9 / queries, sum);
    dbgps_destroy(g);
    return 0;
}
