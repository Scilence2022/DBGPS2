#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Track live allocations so rejected candidates and repeated queries cannot
   leak unnoticed on platforms whose AddressSanitizer lacks leak detection. */
static size_t live_allocations;
static void *tracked_malloc(size_t n)
{
    void *p = malloc(n);
    if (p) ++live_allocations;
    return p;
}
static void *tracked_calloc(size_t n, size_t size)
{
    void *p = calloc(n, size);
    if (p) ++live_allocations;
    return p;
}
static void *tracked_realloc(void *old, size_t n)
{
    int had_pointer = old != NULL;
    void *p = realloc(old, n);
    if (p && !had_pointer) ++live_allocations;
    return p;
}
static void tracked_free(void *p)
{
    if (p) {
        assert(live_allocations > 0);
        --live_allocations;
    }
    free(p);
}

#define malloc tracked_malloc
#define calloc tracked_calloc
#define realloc tracked_realloc
#define free tracked_free
#define main tool_main
#ifdef TEST_DBGPS3
#include "dbgps3_core.c"
#define KC_BITS DBGPS3_KC_BITS
#define canonical dbgps3_min_hash_key
#else
#include "DBGPS2.c"
#define canonical min_hash_key
#endif
#undef main

static const int test_k = 31;
static const uint64_t test_mask = (1ULL << 62) - 1;

static void insert_node(kc_c4x_t *h, uint64_t kmer)
{
    uint64_t hash = hash64(canonical(kmer, test_k), test_mask);
    unsigned partition = hash & ((1U << KC_BITS) - 1);
    int absent;
    khint_t pos = kc_c4_put(h->h[partition], hash & ~((1ULL << KC_BITS) - 1), &absent);
    kh_key(h->h[partition], pos) |= 10;
}

static void add_path(kc_c4x_t *h, uint64_t root, const char *suffix)
{
    insert_node(h, root);
    for (; *suffix; ++suffix) {
        const char *base = strchr("ACGT", *suffix);
        assert(base);
        root = ((root << 2) | (unsigned)(base - "ACGT")) & test_mask;
        insert_node(h, root);
    }
}

static void destroy_table(kc_c4x_t *h)
{
    for (unsigned i = 0; i < 1U << KC_BITS; ++i) kc_c4_destroy(h->h[i]);
    free(h->h);
    free(h);
    assert(live_allocations == 0);
}

#ifdef TEST_DBGPS3
static int query(kc_c4x_t *h, uint64_t root, int path_len, int skip, double ratio,
                 char *out, size_t out_size)
{
    dbgps3_graph_t g = {0};
    g._h = h;
    g.k = test_k;
    g.mask = test_mask;
    return dbgps3_find_path(&g, root, 1, ratio, skip, path_len, 0, out, (int)out_size);
}

static void test_unique_and_ambiguous(void)
{
    kc_c4x_t *h = c4x_init(KC_BITS);
    uint64_t root = dbgps3_actgkmer_uint64("CCTGCAGAGTAGCATGTCACGATCTGACGTAC", test_k);
    /* A splits into two dead ends before C is visited. Writing the next
       generation over the current one used to discard C's surviving path. */
    add_path(h, root, "AA");
    add_path(h, root, "AC");
    add_path(h, root, "CGTACGTAC");
    size_t baseline = live_allocations;
    char out[64];
    for (int i = 0; i < 100; ++i) {
        strcpy(out, "unchanged");
        assert(query(h, root, 8, 0, 2.55, out, sizeof(out)) == 0);
        assert(strcmp(out, "CGTACGTA") == 0);
        assert(live_allocations == baseline);
    }
    add_path(h, root, "CTGACGTAC");
    baseline = live_allocations;
    strcpy(out, "unchanged");
    assert(query(h, root, 8, 0, 2.55, out, sizeof(out)) == -1);
    assert(strcmp(out, "unchanged") == 0);
    assert(live_allocations == baseline);
    assert(query(h, root, -4, 0, 2.55, out, sizeof(out)) == -1);
    assert(query(h, root, 7, 0, 2.55, out, sizeof(out)) == -1);
    assert(query(h, root, 304, 0, 2.55, out, sizeof(out)) == -1);
    assert(query(h, root, 8, 0, 0.0, out, sizeof(out)) == -1);
    assert(query(h, root, 8, 0, NAN, out, sizeof(out)) == -1);
    assert(query(h, root, 8, 0, 2.55, out, 8) == -1);
    assert(strcmp(out, "unchanged") == 0);
    destroy_table(h);
}
#else
static void encoded_strand(unsigned char payload, unsigned char *encoded, char *suffix)
{
    unsigned char raw[4] = {73, payload, 0, 0};
    unsigned short crc = crc16xmodem(raw, 2);
    raw[2] = crc >> 8;
    raw[3] = crc & 255;
    for (int i = 0; i < 4; ++i) encoded[i] = two_to_one_table[raw[i]];
    for (int i = 0; i < 12; ++i)
        suffix[i] = "ACGT"[(encoded[1 + i / 4] >> (6 - 2 * (i % 4))) & 3];
    suffix[12] = 'A'; /* one-base lookahead required by path reconstruction */
    suffix[13] = '\0';
}

static void test_unique_and_ambiguous(void)
{
    kc_c4x_t *h = c4x_init(KC_BITS);
    unsigned char encoded[4];
    char suffix[14];
    encoded_strand(100, encoded, suffix);
    uint64_t root = actgkmer_uint64((unsigned char *)"CCTGCAGAGTAGCATGTCACGATCTGACGTAC", test_k);
    root = (root & ~255ULL) | encoded[0];
    add_path(h, root, "AA");
    add_path(h, root, "AC");
    add_path(h, root, suffix);
    size_t baseline = live_allocations;
    for (int i = 0; i < 100; ++i) {
        unsigned char *paths[Max_Path_Num] = {0};
        int pos = find_paths(h, paths, root, test_mask, 1, 12, test_k, 1, 2.55, 0);
        assert(pos >= 0);
        assert(memcmp(paths[pos], encoded, sizeof(encoded)) == 0);
        free(paths[pos]);
        assert(live_allocations == baseline);
    }
    encoded_strand(101, encoded, suffix);
    add_path(h, root, suffix);
    baseline = live_allocations;
    unsigned char *paths[Max_Path_Num] = {0};
    assert(find_paths(h, paths, root, test_mask, 1, 12, test_k, 1, 2.55, 0) == -1);
    for (int i = 0; i < Max_Path_Num; ++i) assert(paths[i] == NULL);
    assert(live_allocations == baseline);
    destroy_table(h);

    h = c4x_init(KC_BITS);
    suffix[11] = suffix[11] == 'A' ? 'C' : 'A'; /* corrupt the CRC */
    add_path(h, root, suffix);
    baseline = live_allocations;
    assert(find_paths(h, paths, root, test_mask, 1, 12, test_k, 1, 2.55, 0) == -1);
    assert(live_allocations == baseline);
    assert(find_paths(h, paths, root, test_mask, 1, -4, test_k, 1, 2.55, 0) == -1);
    assert(find_paths(h, paths, root, test_mask, 1, 7, test_k, 1, 2.55, 0) == -1);
    assert(find_paths(h, paths, root, test_mask, 1, 300, test_k, 1, 2.55, 0) == -1);
    assert(find_paths(h, paths, root, test_mask, 1, 12, test_k, 1, 0.0, 0) == -1);
    assert(find_paths(h, paths, root, test_mask, 1, 12, test_k, 1, NAN, 0) == -1);
    destroy_table(h);
}

static void custom_strand(const dbgps_codec_t *map, unsigned char payload,
                          unsigned char *packed, char *suffix)
{
    unsigned char raw[12] = {1, 35, 69, 103, 137, payload, 67, 215, 128, 11, 0, 0};
    unsigned short crc = crc16xmodem(raw, 10);
    raw[10] = crc >> 8;
    raw[11] = crc & 255;
    memset(packed, 0, 15);
    /* 12 bytes + four zero bits -> twenty 5-bit words -> sixty bases. */
    for (unsigned block = 0; block < 20; ++block) {
        unsigned value = 0;
        for (unsigned bit = 0; bit < 5; ++bit) {
            unsigned pos = block * 5 + bit;
            value = (value << 1) | (pos < 96 ? ((raw[pos / 8] >> (7 - pos % 8)) & 1) : 0);
        }
        unsigned word = map->encode[value];
        for (unsigned j = 0; j < 3; ++j) {
            unsigned pos = block * 3 + j;
            packed[pos / 4] |= ((word >> (4 - 2 * j)) & 3) << (6 - 2 * (pos % 4));
        }
    }
    for (int pos = 24; pos < 60; ++pos)
        suffix[pos - 24] = "ACGT"[(packed[pos / 4] >> (6 - 2 * (pos % 4))) & 3];
    suffix[36] = 'C';
    suffix[37] = '\0';
}

static void test_custom_ambiguity(void)
{
    dbgps_codec_t map = {0};
    map.bases = 3;
    map.bits = 5;
    map.encode = malloc(32 * sizeof(*map.encode));
    map.decode = malloc(64 * sizeof(*map.decode));
    assert(map.encode && map.decode);
    for (unsigned i = 0; i < 64; ++i) map.decode[i] = UINT32_MAX;
    for (unsigned i = 0; i < 32; ++i) {
        map.encode[i] = (i * 5 + 11) & 63;
        map.decode[map.encode[i]] = i;
    }
    kc_c4x_t *h = c4x_init(KC_BITS);
    unsigned char packed[15];
    char suffix[38];
    uint64_t index;
    assert(dbgps_codec_encode_index(&map, 0x0123456789ULL, 24, &index) == 0);
    uint64_t root = actgkmer_uint64((unsigned char *)"CCTGCAGAGTAGCATGTCACGATCTGACGTAC", test_k);
    root = (root & ~((1ULL << 48) - 1)) | index;
    custom_strand(&map, 100, packed, suffix);
    add_path(h, root, suffix);
    size_t baseline = live_allocations;
    unsigned char *paths[Max_Path_Num] = {0};
    for (int repeat = 0; repeat < 10; ++repeat) {
        int pos = find_paths_mapped(h, paths, root, test_mask, 1, 36, test_k, 24, 2.55, 0, &map);
        assert(pos >= 0);
        assert(memcmp(paths[pos], packed, sizeof(packed)) == 0);
        free(paths[pos]);
        paths[pos] = NULL;
        assert(live_allocations == baseline);
    }
    custom_strand(&map, 101, packed, suffix);
    add_path(h, root, suffix);
    baseline = live_allocations;
    assert(find_paths_mapped(h, paths, root, test_mask, 1, 36, test_k, 24, 2.55, 0, &map) == -1);
    for (int i = 0; i < Max_Path_Num; ++i) assert(paths[i] == NULL);
    assert(live_allocations == baseline);
    /* These allocations are tracked in this translation unit. */
    free(map.encode);
    free(map.decode);
    destroy_table(h);
}
#endif

static void test_branch_limit(void)
{
    kc_c4x_t *h = c4x_init(KC_BITS);
    uint64_t root = 0x123456789abcdULL;
    /* All 4^6 suffixes force more than 1050 simultaneous paths. */
    for (unsigned suffix = 0; suffix < 4096; ++suffix) {
        uint64_t node = root;
        insert_node(h, node);
        for (int shift = 10; shift >= 0; shift -= 2) {
            node = ((node << 2) | ((suffix >> shift) & 3)) & test_mask;
            insert_node(h, node);
        }
    }
    size_t baseline = live_allocations;
    for (int skip = 0; skip <= 8; skip += 8) {
#ifdef TEST_DBGPS3
        char out[64] = "unchanged";
        assert(query(h, root, 8, skip, 2.55, out, sizeof(out)) == -1);
        assert(strcmp(out, "unchanged") == 0);
#else
        unsigned char *paths[Max_Path_Num] = {0};
        assert(find_paths(h, paths, root, test_mask, 1, 8, test_k, 1, 2.55, skip) == -1);
        for (int i = 0; i < Max_Path_Num; ++i) assert(paths[i] == NULL);
#endif
        assert(live_allocations == baseline);
    }
    destroy_table(h);
}

int main(void)
{
    test_unique_and_ambiguous();
#ifndef TEST_DBGPS3
    test_custom_ambiguity();
#endif
    test_branch_limit();
    puts("paths: unique, ambiguous, branching, invalid inputs and allocation checks passed");
    return 0;
}
