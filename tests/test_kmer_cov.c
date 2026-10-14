#include <assert.h>
#include <stdio.h>

/* Compile once per tool to exercise each copy of the coverage lookup. */
#define main tool_main
#include TEST_SOURCE
#undef main

#ifdef TEST_DBGPS3
#define KC_BITS DBGPS3_KC_BITS
#define query_cov kmer_cov_internal
#else
#define query_cov kmer_cov
#endif

int main(void)
{
    const uint64_t mask = (1ULL << 62) - 1;
    kc_c4x_t *h = c4x_init(KC_BITS);
    const unsigned partition = 7;

    /* A 32-slot table has exactly one bitmap word. Probing kh_end here
       previously read beyond that word, unlike some smaller tables. */
    for (unsigned i = 1; i <= 24; ++i) {
        int absent;
        uint64_t key = (uint64_t)i << KC_BITS;
        khint_t pos = kc_c4_put(h->h[partition], key, &absent);
        assert(absent == 1);
        kh_key(h->h[partition], pos) |= 9;
    }
    assert(kh_end(h->h[partition]) == 32);
    assert(kc_c4_get(h->h[partition], (uint64_t)25 << KC_BITS) == kh_end(h->h[partition]));

    for (unsigned i = 1; i <= 128; ++i) {
        uint64_t key = ((uint64_t)i << KC_BITS) | partition;
        assert(query_cov(hash64i(key, mask), mask, h) == (i <= 24 ? 9 : 0));
    }
    assert(query_cov(hash64i(8, mask), mask, h) == 0); /* empty partition */

    /* Make an unchecked kh_end probe fail deterministically even without
       ASan: guard storage is allocated, but lies outside the table capacity. */
    kc_c4_t *table = h->h[partition];
    table->used = realloc(table->used, 2 * sizeof(*table->used));
    table->keys = realloc(table->keys, 33 * sizeof(*table->keys));
    assert(table->used && table->keys);
    table->used[1] = 1;
    kh_key(table, 32) = 123;
    assert(query_cov(hash64i(((uint64_t)25 << KC_BITS) | partition, mask), mask, h) == 0);

    for (unsigned i = 0; i < 1U << KC_BITS; ++i) kc_c4_destroy(h->h[i]);
    free(h->h);
    free(h);
    puts("coverage: present, absent, empty partitions and end-sentinel checks passed");
    return 0;
}
