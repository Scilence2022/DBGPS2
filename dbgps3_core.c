#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "khashl.h"
#include "kthread.h"
#include "dbgps3_core.h"

/* ---- Internal hash table types ---- */

#define kc_c4_eq(a, b) ((a)>>DBGPS3_KC_BITS == (b)>>DBGPS3_KC_BITS)
#define kc_c4_hash(a) ((a)>>DBGPS3_KC_BITS)

KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

/* ---- Translation tables (extern-declared in header) ---- */

const unsigned char dbgps3_seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const unsigned char dbgps3_nt4_seq_table[4] = { 'A', 'C', 'G', 'T' };

/* ---- Internal hash functions ---- */

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

static inline uint64_t hash64i(uint64_t key, uint64_t mask)
{
    uint64_t tmp;
    tmp = (key - (key << 31));    key = (key - (tmp << 31)) & mask;
    tmp = key ^ key >> 28;        key = key ^ tmp >> 28;
    key = (key * 14933078535860113213ull) & mask;
    tmp = key ^ key >> 14; tmp = key ^ tmp >> 14; tmp = key ^ tmp >> 14; key = key ^ tmp >> 14;
    key = (key * 15244667743933553977ull) & mask;
    tmp = key ^ key >> 24;        key = key ^ tmp >> 24;
    tmp = ~key; tmp = ~(key - (tmp << 21)); tmp = ~(key - (tmp << 21)); key = ~(key - (tmp << 21)) & mask;
    return key;
}

/* ---- Partitioned hash table ---- */

typedef struct {
    int p;
    kc_c4_t **h;
} kc_c4x_t;

static kc_c4x_t *c4x_init(int p)
{
    int i;
    kc_c4x_t *h;
    CALLOC(h, 1);
    MALLOC(h->h, 1<<p);
    h->p = p;
    for (i = 0; i < 1<<p; ++i)
        h->h[i] = kc_c4_init();
    return h;
}

static void c4x_destroy(kc_c4x_t *h)
{
    int i;
    if (!h) return;
    for (i = 0; i < 1<<h->p; ++i) kc_c4_destroy(h->h[i]);
    free(h->h);
    free(h);
}

typedef struct {
    int n, m;
    uint64_t *a;
} buf_c4_t;

static inline void c4x_insert_buf(buf_c4_t *buf, int p, uint64_t y)
{
    int pre = y & ((1<<p) - 1);
    buf_c4_t *b = &buf[pre];
    if (b->n == b->m) {
        b->m = b->m < 8? 8 : b->m + (b->m>>1);
        REALLOC(b->a, b->m);
    }
    b->a[b->n++] = y;
}

static void count_seq_buf(buf_c4_t *buf, int k, int p, int len, const char *seq)
{
    int i, l;
    uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = dbgps3_seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) {
            x[0] = (x[0] << 2 | c) & mask;
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
            if (++l >= k) {
                uint64_t y = x[0] < x[1]? x[0] : x[1];
                c4x_insert_buf(buf, p, hash64(y, mask));
            }
        } else l = 0, x[0] = x[1] = 0;
    }
}

/* ---- Pipeline data structures ---- */

typedef struct {
    int k, block_len, n_thread, read_len;
    kseq_t *ks;
    kc_c4x_t *h;
} pldat_t;

typedef struct {
    pldat_t *p;
    int n, m, sum_len, nk;
    int *len;
    char **seq;
    buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid)
{
    stepdat_t *s = (stepdat_t*)data;
    buf_c4_t *b = &s->buf[i];
    kc_c4_t *h = s->p->h->h[i];
    int j, p = s->p->h->p;
    for (j = 0; j < b->n; ++j) {
        khint_t k;
        int absent;
        k = kc_c4_put(h, b->a[j]>>p<<DBGPS3_KC_BITS, &absent);
        if ((kh_key(h, k)&DBGPS3_KC_MAX) < DBGPS3_KC_MAX) ++kh_key(h, k);
    }
}

static void *worker_pipeline(void *data, int step, void *in)
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) {
        int ret;
        stepdat_t *s;
        CALLOC(s, 1);
        s->p = p;
        while ((ret = kseq_read(p->ks)) >= 0) {
            int l = p->ks->seq.l;
            if (l < p->k) continue;
            if (s->n == s->m) {
                s->m = s->m < 16? 16 : s->m + (s->n>>1);
                REALLOC(s->len, s->m);
                REALLOC(s->seq, s->m);
            }
            MALLOC(s->seq[s->n], l);
            memcpy(s->seq[s->n], p->ks->seq.s, l);
            s->len[s->n++] = l;
            s->sum_len += l;
            s->nk += l - p->k + 1;
            if (s->sum_len >= p->block_len)
                break;
        }
        if (s->sum_len == 0) free(s);
        else return s;
    } else if (step == 1) {
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p, m;
        CALLOC(s->buf, n);
        m = (int)(s->nk * 1.2 / n) + 1;
        for (i = 0; i < n; ++i) {
            s->buf[i].m = m;
            MALLOC(s->buf[i].a, m);
        }
        for (i = 0; i < s->n; ++i) {
            if(p->read_len > 0 && p->read_len < s->len[i]){
                count_seq_buf(s->buf, p->k, p->h->p, p->read_len, s->seq[i]);
            }else{
                count_seq_buf(s->buf, p->k, p->h->p, s->len[i], s->seq[i]);
            }
            free(s->seq[i]);
        }
        free(s->seq); free(s->len);
        return s;
    } else if (step == 2) {
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p;
        kt_for(p->n_thread, worker_for, s, n);
        for (i = 0; i < n; ++i) free(s->buf[i].a);
        free(s->buf); free(s);
    }
    return 0;
}

static kc_c4x_t *count_file(const char *fn, int k, int p, int block_size, int n_thread, int read_len)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.read_len = read_len > k ? read_len : k;
    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

static kc_c4x_t *count_file2(const char *fn, kc_c4x_t *hh, int k, int p, int block_size, int n_thread, int read_len)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.read_len = read_len > k ? read_len : k;
    pl.n_thread = n_thread;
    pl.h = hh;
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

/* ---- Internal k-mer coverage query ---- */

static unsigned short kmer_cov_internal(uint64_t kmer, uint64_t mask, const kc_c4x_t *h)
{
    int j, x;
    unsigned short cov = 0;
    uint64_t hash_key = hash64(kmer, mask);
    j = hash_key & ((1<<DBGPS3_KC_BITS) - 1);
    if (kh_size(h->h[j]) < 1) return 0;
    hash_key = hash_key >> DBGPS3_KC_BITS << DBGPS3_KC_BITS;
    x = kc_c4_get(h->h[j], hash_key);
    if (kh_exist(h->h[j], x)) {
        cov = kh_key(h->h[j], x) & DBGPS3_KC_MAX;
    }
    return cov;
}

/* ---- Internal encoding tables ---- */

static const unsigned char two_to_one_table[256] = {
    51, 50, 49, 48, 60, 61, 62, 63, 58, 56, 59, 57, 53, 55, 52, 54,
    35, 34, 33, 32, 44, 45, 46, 47, 42, 40, 43, 41, 37, 39, 36, 38,
    19, 18, 17, 16, 28, 29, 30, 31, 26, 24, 27, 25, 21, 23, 20, 22,
     3,  2,  1,  0, 12, 13, 14, 15, 10,  8, 11,  9,  5,  7,  4,  6,
   195,194,193,192,204,205,206,207,202,200,203,201,197,199,196,198,
   211,210,209,208,220,221,222,223,218,216,219,217,213,215,212,214,
   227,226,225,224,236,237,238,239,234,232,235,233,229,231,228,230,
   243,242,241,240,252,253,254,255,250,248,251,249,245,247,244,246,
   163,162,161,160,172,173,174,175,170,168,171,169,165,167,164,166,
   131,130,129,128,140,141,142,143,138,136,139,137,133,135,132,134,
   179,178,177,176,188,189,190,191,186,184,187,185,181,183,180,182,
   147,146,145,144,156,157,158,159,154,152,155,153,149,151,148,150,
    83, 82, 81, 80, 92, 93, 94, 95, 90, 88, 91, 89, 85, 87, 84, 86,
   115,114,113,112,124,125,126,127,122,120,123,121,117,119,116,118,
    67, 66, 65, 64, 76, 77, 78, 79, 74, 72, 75, 73, 69, 71, 68, 70,
    99, 98, 97, 96,108,109,110,111,106,104,107,105,101,103,100,102
};

/* ---- Internal path node type ---- */

typedef struct path_node_s {
    struct path_node_s *pre_node;
    uint64_t kmer;
    unsigned int cov;
} path_node_t;

static const uint64_t BI_MASK = 3ULL;
static const uint64_t BYT_MASK = 255ULL;

static unsigned char four_kmers_to_byte(uint64_t k1, uint64_t k2, uint64_t k3, uint64_t k4)
{
    unsigned char m;
    m = (unsigned char)(((k1 & BI_MASK)<<6) + ((k2 & BI_MASK)<<4) + ((k3 & BI_MASK)<<2) + ((k4 & BI_MASK)));
    return m;
}

static unsigned char* nodes_path(path_node_t *term_node, unsigned char *path_seq,
                                 int path_len, int km_len, int index_byte_num)
{
    path_node_t *p_node;
    p_node = term_node;
    int path_seq_p = (path_len>>2) - 1 + index_byte_num, nd_p = 0;
    uint64_t kms[4];
    kms[nd_p] = p_node->kmer;

    while(p_node->pre_node && path_seq_p >= index_byte_num){
        p_node = p_node->pre_node;
        kms[nd_p] = p_node->kmer;

        if(nd_p >= 3){
            *(path_seq + path_seq_p) = four_kmers_to_byte(kms[3], kms[2], kms[1], kms[0]);
            nd_p = 0;
            path_seq_p--;
        } else { nd_p++; }
    }

    while(p_node->pre_node){
        p_node = p_node->pre_node;
    }
    uint64_t index_kmer = p_node->kmer;
    int i;
    for(i=index_byte_num-1; i>=0; i--){
        *(path_seq + i) = (unsigned char) (index_kmer & BYT_MASK);
        index_kmer = index_kmer>>8;
    }

    return path_seq;
}

static void intseq_actg(const unsigned char *int_seq, unsigned char *actg_seq, int8_t int_seq_len)
{
    int p = 0;
    unsigned long n, mask;
    unsigned char n1, n2, n3, n4;

    while(p < int_seq_len){
        unsigned long tmp = 0;
        int i2;
        /* Decode one byte to 4 bases */
        unsigned char bv = int_seq[p];
        for(i2=0; i2<4; i2++){
            int m = bv & 3;
            bv = bv >> 2;
            tmp = tmp >> 8;
            tmp = tmp + ((unsigned long)dbgps3_nt4_seq_table[m] << 24);
        }
        mask = (1<<8) - 1;
        n = tmp;
        n4 = (n & mask);
        mask = mask << 8;
        n3 = (n & mask)>>8;
        mask = mask << 8;
        n2 = (n & mask)>>16;
        mask = mask << 8;
        n1 = (n & mask)>>24;

        actg_seq[p<<2] = n1;
        actg_seq[(p<<2) + 1] = n2;
        actg_seq[(p<<2) + 2] = n3;
        actg_seq[(p<<2) + 3] = n4;
        p = p + 1;
    }
}

/* ---- Graph struct definition ---- */
/* The full definition is in dbgps3_core.h; here we just use it directly. */
/* The internal hash table is accessed via (kc_c4x_t*)g->_h */

/* ============================================================ */
/* Public API implementations                                    */
/* ============================================================ */

dbgps3_graph_t *dbgps3_load(int k, int p, int n_thread, int block_size,
                             int read_len, int n_files, char **files)
{
    if (n_files < 1 || !files) return NULL;

    dbgps3_graph_t *g;
    CALLOC(g, 1);
    g->k = k;
    g->p = p;
    g->mask = (1ULL<<k*2) - 1;
    g->n_thread = n_thread;
    g->block_size = block_size;
    g->read_len = read_len > k ? read_len : k;

    /* Default index/payload parameters — can be overridden by caller */
    g->index_len = 16;
    g->index_byte_num = g->index_len >> 2;
    g->path_len = 140; /* default payload_len */

    fprintf(stderr, "[dbgps3_core] Counting k-mers from file 1 ......\n");
    g->_h = count_file(files[0], k, p, block_size, n_thread, read_len);
    if (!g->_h) {
        free(g);
        return NULL;
    }

    int i;
    for (i = 1; i < n_files; i++) {
        fprintf(stderr, "[dbgps3_core] Counting k-mers from file %d ......\n", i+1);
        g->_h = count_file2(files[i], (kc_c4x_t*)g->_h, k, p, block_size, n_thread, read_len);
    }

    return g;
}

int dbgps3_load_more(dbgps3_graph_t *g, int n_files, char **files)
{
    if (!g || n_files < 1 || !files) return -1;
    int i;
    for (i = 0; i < n_files; i++) {
        fprintf(stderr, "[dbgps3_core] Loading additional file %d ......\n", i+1);
        g->_h = count_file2(files[i], (kc_c4x_t*)g->_h, g->k, g->p, g->block_size, g->n_thread, g->read_len);
        if (!g->_h) return -1;
    }
    return 0;
}

void dbgps3_destroy(dbgps3_graph_t *g)
{
    if (!g) return;
    c4x_destroy((kc_c4x_t*)g->_h);
    free(g);
}

unsigned short dbgps3_kmer_cov(const dbgps3_graph_t *g, uint64_t kmer)
{
    if (!g) return 0;
    return kmer_cov_internal(kmer, g->mask, (const kc_c4x_t*)g->_h);
}

int dbgps3_find_path(const dbgps3_graph_t *g, uint64_t start_kmer,
                     int cut_cov, double max_cov_ratio, int skip_ratios,
                     int path_len, int index_byte_num,
                     char *result_buf, int result_buf_len)
{
    if (!g || !result_buf) return -1;

    int k = g->k;
    uint64_t mask = g->mask;
    const kc_c4x_t *h = (const kc_c4x_t*)g->_h;
    int seq_len = (path_len>>2) + index_byte_num;

    /* Check result buffer is large enough */
    int needed = path_len + index_byte_num * 4 + 1;
    if (result_buf_len < needed) return -1;

    /* Allocate workspace locally — no global variables, thread-safe */
    path_node_t *all_nodes = (path_node_t *)malloc(
        DBGPS3_MAX_PATH_NUM * DBGPS3_MAX_PATH_LEN * sizeof(path_node_t));
    path_node_t **term_nodes = (path_node_t **)malloc(
        (DBGPS3_MAX_PATH_NUM + 1) * sizeof(path_node_t *));
    if (!all_nodes || !term_nodes) {
        free(all_nodes);
        free(term_nodes);
        return -1;
    }

    all_nodes[0].kmer = start_kmer;
    all_nodes[0].cov = kmer_cov_internal(dbgps3_min_hash_key(start_kmer, k), mask, h);
    int a_cov = all_nodes[0].cov / max_cov_ratio;
    int b_cov = all_nodes[0].cov * max_cov_ratio;
    if (a_cov < cut_cov + 1) a_cov = cut_cov + 1;

    term_nodes[0] = &all_nodes[0];
    int term_node_num = 1, all_node_num = 1, p = 0;
    path_node_t *p_node;
    path_node_t next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0;
    next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    int new_term_node_num;

    while(p <= path_len && term_node_num > 0 && term_node_num <= DBGPS3_MAX_PATH_NUM){
        int i;
        new_term_node_num = 0;
        int n_p = 0;
        while(n_p < term_node_num && new_term_node_num <= DBGPS3_MAX_PATH_NUM){
            p_node = term_nodes[n_p];
            a_cov = p_node->cov / max_cov_ratio;
            b_cov = p_node->cov * max_cov_ratio;
            if(a_cov < cut_cov + 1) a_cov = cut_cov + 1;

            next_nodes[0].kmer = ((*p_node).kmer<<2) & mask;
            next_nodes[1].kmer = (((*p_node).kmer<<2) + 1) & mask;
            next_nodes[2].kmer = (((*p_node).kmer<<2) + 2) & mask;
            next_nodes[3].kmer = (((*p_node).kmer<<2) + 3) & mask;

            next_nodes[0].cov = kmer_cov_internal(dbgps3_min_hash_key(next_nodes[0].kmer, k), mask, h);
            next_nodes[1].cov = kmer_cov_internal(dbgps3_min_hash_key(next_nodes[1].kmer, k), mask, h);
            next_nodes[2].cov = kmer_cov_internal(dbgps3_min_hash_key(next_nodes[2].kmer, k), mask, h);
            next_nodes[3].cov = kmer_cov_internal(dbgps3_min_hash_key(next_nodes[3].kmer, k), mask, h);

            for(i=0; i<4; i++){
                if (p < skip_ratios) {
                    if(next_nodes[i].cov >= cut_cov + 1) {
                        all_nodes[all_node_num].pre_node = p_node;
                        all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                        all_nodes[all_node_num].cov = next_nodes[i].cov;
                        term_nodes[new_term_node_num] = &all_nodes[all_node_num];
                        new_term_node_num++;
                        all_node_num++;
                    }
                } else {
                    if(next_nodes[i].cov >= a_cov && next_nodes[i].cov <= b_cov) {
                        all_nodes[all_node_num].pre_node = p_node;
                        all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                        all_nodes[all_node_num].cov = next_nodes[i].cov;
                        term_nodes[new_term_node_num] = &all_nodes[all_node_num];
                        new_term_node_num++;
                        all_node_num++;
                    }
                }
            }
            n_p++;
        }
        term_node_num = new_term_node_num;
        ++p;
    }

    int ret;
    if(term_node_num == 1){
        unsigned char *path_buf;
        CALLOC(path_buf, seq_len);
        nodes_path(term_nodes[0], path_buf, path_len, k, index_byte_num);

        /* Decode internal seq to ACGT */
        unsigned char *dna_buf = (unsigned char *)result_buf;
        intseq_actg(path_buf, dna_buf, (path_len>>2) + index_byte_num);
        dna_buf[path_len + index_byte_num*4] = '\0';

        free(path_buf);
        ret = 0;
    } else {
        ret = -1;
    }

    free(all_nodes);
    free(term_nodes);
    return ret;
}

/* ---- Utility function implementations ---- */

uint64_t dbgps3_comp_rev(uint64_t x, unsigned char km_len)
{
    uint64_t y = 0;
    int i, c;
    for(i=0; i<km_len; i++){
        y = y << 2;
        c = x & 3ULL;
        y = y | (3-c);
        x = x >> 2;
    }
    return y;
}

uint64_t dbgps3_min_hash_key(uint64_t x, unsigned char km_len)
{
    uint64_t y = dbgps3_comp_rev(x, km_len);
    return (x < y) ? x : y;
}

uint64_t dbgps3_actgkmer_uint64(const char *seq, int km_len)
{
    int i;
    uint64_t x, mask = (1ULL<<km_len*2) - 1;
    for (i = 0, x = 0; i < km_len; ++i) {
        int c = dbgps3_seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) {
            x = (x << 2 | c) & mask;
        } else x = 0;
    }
    return x;
}

uint64_t dbgps3_actgkmer_hashkey(const char *seq, int km_len)
{
    int i;
    uint64_t x[2], y, shift = (km_len - 1) * 2, mask = (1ULL<<km_len*2) - 1;
    for (i = 0, x[0] = x[1] = 0; i < km_len; ++i) {
        int c = dbgps3_seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) {
            x[0] = (x[0] << 2 | c) & mask;
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
        } else x[0] = x[1] = 0;
    }
    y = x[0] < x[1]? x[0] : x[1];
    return y;
}

uint64_t dbgps3_primer_uint64(const char *pr, int pr_len, int index_len, int km_len)
{
    unsigned char *seq;
    uint64_t pr_int64;
    int i = pr_len + index_len - km_len, j = 0;

    MALLOC(seq, km_len);

    while(i < pr_len){
        *(seq + j) = *(pr + i);
        i++;
        j++;
    }

    while(j < km_len){
        *(seq + j) = 'A';
        j++;
    }

    pr_int64 = dbgps3_actgkmer_uint64((const char *)seq, km_len);
    free(seq);
    return pr_int64;
}

uint64_t dbgps3_index_encode(uint64_t index, unsigned int index_len)
{
    uint64_t mask = 255ULL;
    unsigned int index_byte_num = index_len>>2;
    unsigned int i;
    uint64_t m, y = 0ULL;

    for(i=0; i < index_byte_num; i++){
        m = index & mask;
        m = (uint64_t)two_to_one_table[m];
        m = m<<((index_byte_num-1)*8);
        y = y>>8;
        y = y | m;
        index = index>>8;
    }
    return y;
}

void dbgps3_uint64_acgt(uint64_t key, char *seq, int km_len)
{
    int p = km_len - 1;
    while (p >= 0){
        int n = key % 4;
        key = key >> 2;
        seq[p] = dbgps3_nt4_seq_table[n];
        p = p - 1;
    }
}

/* ---- Graph accessors ---- */

int dbgps3_get_k(const dbgps3_graph_t *g) { return g ? g->k : 0; }
uint64_t dbgps3_get_mask(const dbgps3_graph_t *g) { return g ? g->mask : 0; }
int dbgps3_get_index_len(const dbgps3_graph_t *g) { return g ? g->index_len : 0; }
int dbgps3_get_index_byte_num(const dbgps3_graph_t *g) { return g ? g->index_byte_num : 0; }
int dbgps3_get_path_len(const dbgps3_graph_t *g) { return g ? g->path_len : 0; }
int dbgps3_get_p(const dbgps3_graph_t *g) { return g ? g->p : 0; }

long dbgps3_get_kmer_count(const dbgps3_graph_t *g)
{
    if (!g || !g->_h) return 0;
    const kc_c4x_t *h = (const kc_c4x_t*)g->_h;
    long count = 0;
    int i;
    for (i = 0; i < 1<<h->p; ++i)
        count += kh_size(h->h[i]);
    return count;
}
