#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static int fail_after = -1;
static size_t live;
static void *test_calloc(size_t count, size_t size)
{
    if (fail_after == 0) return NULL;
    if (fail_after > 0) --fail_after;
    void *p = calloc(count, size);
    if (p) ++live;
    return p;
}
static void test_free(void *p)
{
    if (p) { assert(live > 0); --live; }
    free(p);
}
#define calloc test_calloc
#define free test_free
#include "dbgps_bloom.h"

static uint64_t checksum(const dbgps_bloom_t *b)
{
    uint64_t sum = 0;
    for (size_t i = 0; i < b->n_bytes; ++i) sum = sum * 131 + b->bits[i];
    return sum;
}

int main(void)
{
    const char *invalid[] = {"", "-1", "+12", " 12", "12 ", "1.2", "65",
                             "12x", "99999999999999999999999999999"};
    int bits = -1;
    for (size_t i = 0; i < sizeof(invalid) / sizeof(*invalid); ++i)
        assert(dbgps_bloom_parse_bits(invalid[i], &bits) == -1);
    assert(dbgps_bloom_parse_bits("0", &bits) == 0 && bits == 0);
    assert(dbgps_bloom_parse_bits("12", &bits) == 0 && bits == 12);
    assert(dbgps_bloom_parse_bits("64", &bits) == 0 && bits == 64);
    assert(!dbgps_bloom_create(10, 0));
    assert(!dbgps_bloom_create(10, 65));
    assert(!dbgps_bloom_create(UINT64_MAX, 64));
    for (int failure = 0; failure < 2; ++failure) {
        fail_after = failure;
        assert(!dbgps_bloom_create(100, 12) && live == 0);
    }
    fail_after = -1;
    dbgps_bloom_t *b = dbgps_bloom_create(0, 12);
    assert(b && b->n_bytes == 64 && !((uintptr_t)b->bits & 63));
    assert(!dbgps_bloom_contains(b, 0));
    dbgps_bloom_insert(b, 0);
    dbgps_bloom_insert(b, UINT64_MAX);
    assert(dbgps_bloom_contains(b, 0) && dbgps_bloom_contains(b, UINT64_MAX));
    dbgps_bloom_destroy(b);

    const unsigned n = 100000;
    b = dbgps_bloom_create(n, 12);
    assert(b && b->n_bytes * 8 >= (size_t)n * 12);
    assert(b->n_bytes * 8 < (size_t)n * 24);
    for (unsigned i = 0; i < n; ++i) dbgps_bloom_insert(b, i);
    uint64_t before = checksum(b);
    for (unsigned i = 0; i < n; ++i) assert(dbgps_bloom_contains(b, i));
    unsigned positives = 0;
    for (unsigned i = n; i < 2 * n; ++i) positives += dbgps_bloom_contains(b, i);
    /* Broad deterministic bound catches broken hash/block selection. */
    assert(positives < n / 50 && checksum(b) == before);
    printf("Bloom: no false negatives; %u/%u absent keys positive; sizing, parsing and allocation failures passed\n",
           positives, n);
    dbgps_bloom_destroy(b);
    dbgps_bloom_destroy(NULL);
    assert(live == 0 && dbgps_bloom_contains(NULL, 1));
    return 0;
}
