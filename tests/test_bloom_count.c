#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define main tool_main
#include TEST_SOURCE
#undef main

#ifdef TEST_DBGPS2_RPC
#define KC_BITS DBGPS_KC_BITS
#define KC_MAX DBGPS_KC_MAX
#define query_cov kmer_cov_internal
#else
#define query_cov kmer_cov
#endif

static uint64_t rng = 91;
static uint64_t next_random(void)
{
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return rng;
}

static void compare_counts(kc_c4x_t *exact, kc_c4x_t *filtered, uint64_t mask)
{
    size_t retained = 0, actual = 0, singletons = 0;
    assert(c4x_set_bloom(filtered, 12) == 0);
    for (int i = 0; i < 1 << exact->p; ++i) {
        kc_c4_t *table = exact->h[i];
        actual += kh_size(filtered->h[i]);
        for (khint_t j = 0; j < kh_end(table); ++j) {
            if (!kh_exist(table, j)) continue;
            unsigned cov = kh_key(table, j) & KC_MAX;
            uint64_t hash = (kh_key(table, j) >> KC_BITS << exact->p) | (uint64_t)i;
            assert(query_cov(hash64i(hash, mask), mask, filtered) == (cov >= 2 ? cov : 0));
            retained += cov >= 2;
            singletons += cov == 1;
        }
    }
    assert(retained == actual);
    printf("counting: %zu retained k-mers counted exactly; %zu singletons removed\n", retained, singletons);
}

#ifdef TEST_DBGPS2_RPC
static void test_load_more(char **files, const char *missing)
{
    FILE *fp = fopen(files[0], "w");
    assert(fp);
    fprintf(fp, ">once\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
                ">twice1\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
                ">twice2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n");
    fclose(fp);
    fp = fopen(files[1], "w");
    assert(fp);
    fprintf(fp, ">reverse\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
    fclose(fp);
    dbgps_graph_t *g = dbgps_load_filtered(31, KC_BITS, 3, 40, 80, 1, files, 23);
    assert(g && dbgps_set_bloom(g, 12) == 0);
    assert(!strcmp(dbgps_redundancy_name(&g->redundancy), "crc16"));
    uint64_t a = dbgps_actgkmer_hashkey("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 31);
    uint64_t c = dbgps_actgkmer_hashkey("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", 31);
    assert(dbgps_kmer_cov(g, a) == 0 && dbgps_kmer_cov(g, c) == 2);
    assert(dbgps_load_more(g, 1, files + 1) == 0);
    assert(!strcmp(dbgps_redundancy_name(&g->redundancy), "crc16"));
    assert(dbgps_kmer_cov(g, a) == 2 && dbgps_kmer_cov(g, c) == 2);
    assert(g->bloom_bits_per_key == 12 && dbgps_bloom_bytes(g) > 0);
    void *before = g->_h;
    char *bad_files[] = {files[1], (char*)missing};
    assert(dbgps_load_more(g, 2, bad_files) == -1 && g->_h == before);
    assert(dbgps_kmer_cov(g, a) == 2 && dbgps_kmer_cov(g, c) == 2);
    dbgps_destroy(g);

    /* Ordinary counting also preserves singleton counts under query filtering.
     * A failed partial append must disable the stale filter and keep g usable. */
    g = dbgps_load(31, KC_BITS, 1, 40, 80, 1, files);
    assert(g && dbgps_set_bloom(g, 12) == 0);
    assert(!strcmp(dbgps_redundancy_name(&g->redundancy), "crc16"));
    assert(dbgps_set_redundancy(g, "none", NULL, 0) == 0);
    assert(dbgps_kmer_cov(g, a) == 1);
    assert(dbgps_load_more(g, 2, bad_files) == -1);
    assert(dbgps_kmer_cov(g, a) == 2 && !dbgps_bloom_bytes(g));
    assert(dbgps_load_more(g, 1, files + 1) == 0);
    assert(!dbgps_redundancy_enabled(&g->redundancy));
    assert(dbgps_kmer_cov(g, a) == 3 && dbgps_bloom_bytes(g));
    dbgps_destroy(g);
    puts("counting API: append restores prior singletons; rebuild and failure ownership passed");
}
#endif

int main(int argc, char **argv)
{
    assert(argc == 2);
    char paths[2][2048], missing[2048];
    snprintf(paths[0], sizeof(paths[0]), "%s/count1.fa", argv[1]);
    snprintf(paths[1], sizeof(paths[1]), "%s/count2.fa", argv[1]);
    snprintf(missing, sizeof(missing), "%s/missing.fa", argv[1]);
    char *files[] = {paths[0], paths[1]};
    FILE *fp[2] = {fopen(files[0], "w"), fopen(files[1], "w")};
    assert(fp[0] && fp[1]);
    for (int n = 0; n < 3000; ++n) {
        char dna[81], reverse[81];
        for (int j = 0; j < 80; ++j) dna[j] = "ACGT"[next_random() & 3];
        if (n % 13 == 0) dna[12] = 'N';
        dna[80] = reverse[80] = 0;
        const char *alphabet = "ACGT";
        for (int j = 0; j < 80; ++j) {
            const char *base = strchr(alphabet, dna[79 - j]);
            reverse[j] = base ? "TGCA"[base - alphabet] : 'N';
        }
        fprintf(fp[0], ">%d\n%s\n", n, dna);
        /* Reverse-complement copies exercise canonicalization across files. */
        if (n % 3 == 0) fprintf(fp[1], ">%d\n%s\n", n, reverse);
    }
    for (int i = 0; i < 20000; ++i)
        fprintf(fp[i & 1], ">saturation%d\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n", i);
    fclose(fp[0]); fclose(fp[1]);
    for (int k = 7; k <= 31; k += 24) {
        uint64_t mask = (UINT64_C(1) << (2 * k)) - 1;
        kc_c4x_t *exact = count_file(files[0], k, KC_BITS, 10000, 3, 80);
        assert(exact && count_file2(files[1], exact, k, KC_BITS, 10000, 3, 80));
        for (int threads = 1; threads <= 3; threads += 2) {
            kc_c4x_t *filtered = count_files_filtered(2, files, k, KC_BITS, 10000, threads, 80, 23);
            assert(filtered);
            compare_counts(exact, filtered, mask);
            assert(query_cov(0, mask, filtered) == KC_MAX);
            c4x_destroy(filtered);
        }
        c4x_destroy(exact);
    }
    /* Force a false-positive singleton candidate; the recount must remove it. */
    kc_c4x_t *h = c4x_init(KC_BITS);
    h->count_pass = 1;
    h->seen = calloc(1 << KC_BITS, sizeof(*h->seen));
    assert(h->seen);
    h->seen[0] = dbgps_bloom_create(1, 12);
    assert(h->seen[0]);
    memset(h->seen[0]->bits, 255, h->seen[0]->n_bytes);
    c4x_count_key(h, 0, 0);
    assert(kh_size(h->h[0]) == 1);
    h->count_pass = 2;
    c4x_count_key(h, 0, 0);
    c4x_prune_singletons(h);
    assert(kh_size(h->h[0]) == 0);
    c4x_destroy(h);
#ifdef TEST_DBGPS2_RPC
    test_load_more(files, missing);
#endif
    puts("counting: two files, thread counts, reverse complements, Ns, saturation and Bloom false positives passed");
    return 0;
}
