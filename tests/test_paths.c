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
#ifdef TEST_DBGPS2_RPC
#include "dbgps_core.c"
#define KC_BITS DBGPS_KC_BITS
#define canonical dbgps_min_hash_key
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

#ifdef TEST_DBGPS2_RPC
static int query(kc_c4x_t *h, uint64_t root, int path_len, int skip, double ratio,
                 char *out, size_t out_size)
{
    dbgps_graph_t g = {0};
    g._h = h;
    g.k = test_k;
    g.mask = test_mask;
    return dbgps_find_path(&g, root, 1, ratio, skip, path_len, 0, out, (int)out_size);
}

static void test_unique_and_ambiguous(void)
{
    kc_c4x_t *h = c4x_init(KC_BITS);
    uint64_t root = dbgps_actgkmer_uint64("CCTGCAGAGTAGCATGTCACGATCTGACGTAC", test_k);
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

static void test_mapped_paths(void)
{
    dbgps_graph_t g = {0};
    g.k = test_k;
    g.mask = test_mask;
    uint64_t index;
    assert(dbgps_index_encode(0x0123456789abULL, 24) == 0x0123456789abULL);
    assert(dbgps_encode_index(&g, 0x0123456789abULL, 24, &index) == 0);
    assert(index == 0x0123456789abULL);
    assert(dbgps_encode_index(&g, 1ULL << 48, 24, &index) == -1);

    g.codec.bases = 3;
    g.codec.bits = 5;
    g.codec.encode = malloc(32 * sizeof(*g.codec.encode));
    g.codec.decode = malloc(64 * sizeof(*g.codec.decode));
    assert(g.codec.encode && g.codec.decode);
    for (unsigned i = 0; i < 64; ++i) g.codec.decode[i] = i < 32 ? i : UINT32_MAX;
    for (unsigned i = 0; i < 32; ++i) g.codec.encode[i] = i;
    assert(dbgps_encode_index(&g, 0x0123456789ULL, 24, &index) == 0);
    uint64_t root = dbgps_actgkmer_uint64("CCTGCAGAGTAGCATGTCACGATCTGACGTAC", test_k);
    root = (root & ~((1ULL << 48) - 1)) | index;
    kc_c4x_t *h = c4x_init(KC_BITS);
    g._h = h;
    char suffix[26];
    for (unsigned block = 0; block < 8; ++block) {
        unsigned value = (block * 7 + 3) % 32;
        for (unsigned pos = 0; pos < 3; ++pos)
            suffix[block * 3 + pos] = "ACGT"[(value >> (4 - 2 * pos)) & 3];
    }
    suffix[24] = 'C';
    suffix[25] = '\0';
    add_path(h, root, suffix);
    char forbidden[sizeof(suffix)];
    strcpy(forbidden, suffix);
    forbidden[0] = 'T'; /* word >= 48, forbidden by this codebook */
    add_path(h, root, forbidden);

    char expected[49], out[64];
    dbgps_uint64_acgt(index, expected, 24);
    memcpy(expected + 24, suffix, 24);
    expected[48] = '\0';
    size_t baseline = live_allocations;
    strcpy(out, "unchanged");
    assert(dbgps_find_path(&g, root, 1, 2.55, 0, 24, 5, out, 48) == -1);
    assert(!strcmp(out, "unchanged"));
    assert(live_allocations == baseline);
    for (int repeat = 0; repeat < 10; ++repeat) {
        strcpy(out, "unchanged");
        assert(dbgps_find_path(&g, root, 1, 2.55, 0, 24, 5, out, sizeof(out)) == 0);
        assert(!strcmp(out, expected));
        assert(live_allocations == baseline);
    }
    suffix[23] = suffix[23] == 'A' ? 'C' : 'A'; /* another allowed final word */
    add_path(h, root, suffix);
    baseline = live_allocations;
    strcpy(out, "unchanged");
    assert(dbgps_find_path(&g, root, 1, 2.55, 0, 24, 5, out, sizeof(out)) == -1);
    assert(!strcmp(out, "unchanged"));
    assert(live_allocations == baseline);
    free(g.codec.encode);
    free(g.codec.decode);
    destroy_table(h);
}
#else
static void encoded_strand(unsigned char payload, unsigned char *encoded, char *suffix)
{
    unsigned char raw[4] = {73, payload, 0, 0};
    unsigned short crc = crc16xmodem(raw, 2);
    raw[2] = crc >> 8;
    raw[3] = crc & 255;
    memcpy(encoded, raw, sizeof(raw)); /* Native A=00,C=01,G=10,T=11 packing. */
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

static void test_builtin_encoding(void)
{
    /* Integer indices must retain their bit patterns, including leading zero
       bytes and bits beyond the old uint32_t index limit. */
    assert(index_encode(0x01234567ULL, 16) == 0x01234567ULL);
    assert(index_encode(0xff01234567ULL, 16) == 0x01234567ULL);
    assert(index_encode(0x0123456789abULL, 24) == 0x0123456789abULL);
    assert(index_encode(UINT64_MAX, 0) == 0);
    assert(index_encode(UINT64_MAX, 32) == UINT64_MAX);
    for (unsigned byte = 0; byte < 256; ++byte) {
        unsigned char dna[4];
        for (unsigned pos = 0; pos < 4; ++pos)
            dna[pos] = "ACGT"[(byte >> (6 - 2 * pos)) & 3];
        assert(actgkmer_uint64(dna, 4) == byte);
    }
}
#endif

static int protected_query(kc_c4x_t *h, uint64_t root, const dbgps_redundancy_t *c,
                           const unsigned char *expected, unsigned length)
{
#ifdef TEST_DBGPS2_RPC
    dbgps_graph_t g = {0};
    g._h = h;
    g.k = test_k;
    g.mask = test_mask;
    g.redundancy = *c;
    char out[160] = "unchanged";
    int rc = dbgps_find_path(&g, root, 1, 2.55, 0, 24, 4, out, sizeof(out));
    if (rc == 0) {
        assert(strlen(out) == length * 4);
        for (unsigned i = 0; i < length * 4; ++i)
            assert(out[i] == "ACGT"[(expected[i / 4] >> (6 - 2 * (i % 4))) & 3]);
    } else assert(!strcmp(out, "unchanged"));
    return rc;
#else
    unsigned char *paths[Max_Path_Num] = {0};
    int pos = find_paths_protected(h, paths, root, test_mask, 1, 24, test_k, 16, 2.55, 0, NULL, c);
    if (pos < 0) return -1;
    assert(!memcmp(paths[pos], expected, length));
    free(paths[pos]);
    return 0;
#endif
}

static void add_encoded_path(kc_c4x_t *h, uint64_t root, const unsigned char *word, unsigned length)
{
    char suffix[160];
    for (unsigned i = 16; i < length * 4; ++i)
        suffix[i - 16] = "ACGT"[(word[i / 4] >> (6 - 2 * (i % 4))) & 3];
    suffix[length * 4 - 16] = 'C'; /* historical primer lookahead */
    suffix[length * 4 - 15] = '\0';
    add_path(h, root, suffix);
}

static void test_protected_paths(void)
{
    const char *names[] = {"none", "crc16", "crc32", "crc32c", "rs:4", "hamming"};
    for (unsigned mode = 0; mode < sizeof(names) / sizeof(*names); ++mode) {
        dbgps_redundancy_t c = {0};
        assert(dbgps_redundancy_init(&c, names[mode], NULL, 0) == 0);
        unsigned char expected[32] = {0x01,0x35,0xab,0x79,0x26,0x93,0x4d,0xb8,0x65,0xf1};
        size_t parity;
        assert(dbgps_redundancy_size(&c, 10, &parity) == 0);
        unsigned length = 10 + parity;
        assert(dbgps_redundancy_encode(&c, expected, 10, expected, length) == 0);
        uint64_t root = 0x19358bcULL << 32 | 0x0135ab79ULL;
        kc_c4x_t *h = c4x_init(KC_BITS);
        unsigned char damaged[32];
        memcpy(damaged, expected, length);
        int repairs = mode >= 4;
        if (mode) damaged[5] ^= 0x80;
        add_encoded_path(h, root, damaged, length);
        size_t baseline = live_allocations;
        int wanted = mode && !repairs ? -1 : 0;
        for (unsigned repeat = 0; repeat < 10; ++repeat) {
            assert(protected_query(h, root, &c, expected, length) == wanted);
            assert(live_allocations == baseline);
        }
        add_encoded_path(h, root, expected, length);
        baseline = live_allocations;
        assert(protected_query(h, root, &c, expected, length) == 0);
        assert(live_allocations == baseline);
        if (repairs) {
            memcpy(damaged, expected, length);
            damaged[length - 1] ^= 0x80; /* parity-only errors also converge */
            add_encoded_path(h, root, damaged, length);
            baseline = live_allocations;
            assert(protected_query(h, root, &c, expected, length) == 0);
            assert(live_allocations == baseline);
        }
        expected[4] ^= 0x20;
        assert(dbgps_redundancy_encode(&c, expected, 10, expected, length) == 0);
        add_encoded_path(h, root, expected, length);
        baseline = live_allocations;
        assert(protected_query(h, root, &c, expected, length) == -1);
        assert(live_allocations == baseline);
        destroy_table(h);
    }
#ifdef TEST_DBGPS2_RPC
    dbgps_graph_t g = {0};
    g.k = 31; g.index_len = 16; g.path_len = 140;
    char error[256];
    assert(dbgps_set_redundancy(&g, "rs:4", error, sizeof(error)) == 0);
    assert(dbgps_set_redundancy(&g, "rs:254", error, sizeof(error)) == -1);
    assert(g.redundancy.parity == 4);
    assert(dbgps_set_redundancy(&g, "bad", error, sizeof(error)) == -1);
    assert(g.redundancy.parity == 4);
    assert(dbgps_set_redundancy(&g, "none", error, sizeof(error)) == 0);
#endif
}

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
#ifdef TEST_DBGPS2_RPC
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
#ifndef TEST_DBGPS2_RPC
    test_builtin_encoding();
    test_custom_ambiguity();
#else
    test_mapped_paths();
#endif
    test_branch_limit();
    test_protected_paths();
    puts("paths: unique, ambiguous, branching, invalid inputs and allocation checks passed");
    return 0;
}
