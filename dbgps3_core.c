#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "khashl.h"
#include "kthread.h"
#include "dbgps3_core.h"
#include "dbgps_bloom.h"

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

#define DBGPS_COUNT_BITS DBGPS3_KC_BITS
#define DBGPS_COUNT_MAX DBGPS3_KC_MAX
#include "dbgps_bloom_graph.h"

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
    int read_error;
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
    for (int j = 0; j < b->n; ++j)
        c4x_count_key(s->p->h, (int)i, b->a[j]);
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
        if (ret < -1) p->read_error = 1;
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
    pldat_t pl = {0};
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.read_len = read_len > k ? read_len : k;
    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    if (!pl.h) { kseq_destroy(pl.ks); gzclose(fp); return NULL; }
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    int gzip_error;
    gzerror(fp, &gzip_error);
    kseq_destroy(pl.ks);
    int close_error = gzclose(fp);
    if (pl.read_error || (gzip_error != Z_OK && gzip_error != Z_STREAM_END) || close_error != Z_OK) {
        c4x_destroy(pl.h);
        return NULL;
    }
    return pl.h;
}

static kc_c4x_t *count_file2(const char *fn, kc_c4x_t *hh, int k, int p, int block_size, int n_thread, int read_len)
{
    pldat_t pl = {0};
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    pl.read_len = read_len > k ? read_len : k;
    pl.n_thread = n_thread;
    pl.h = hh;
    c4x_set_bloom(pl.h, 0);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    int gzip_error;
    gzerror(fp, &gzip_error);
    kseq_destroy(pl.ks);
    int close_error = gzclose(fp);
    if (pl.read_error || (gzip_error != Z_OK && gzip_error != Z_STREAM_END) || close_error != Z_OK) {
        return NULL; /* caller owns the possibly partially updated table */
    }
    return pl.h;
}

/* ---- Internal k-mer coverage query ---- */

static unsigned short kmer_cov_internal(uint64_t kmer, uint64_t mask, const kc_c4x_t *h)
{
    int j, x;
    unsigned short cov = 0;
    uint64_t hash_key = hash64(kmer, mask);
    if (h->bloom && !dbgps_bloom_contains(h->bloom, hash_key)) return 0;
    j = hash_key & ((1<<DBGPS3_KC_BITS) - 1);
    if (kh_size(h->h[j]) < 1) return 0;
    hash_key = hash_key >> DBGPS3_KC_BITS << DBGPS3_KC_BITS;
    x = kc_c4_get(h->h[j], hash_key);
    if (x != kh_end(h->h[j]) && kh_exist(h->h[j], x)) {
        cov = kh_key(h->h[j], x) & DBGPS3_KC_MAX;
    }
    return cov;
}

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

static void nodes_path_packed(path_node_t *node, unsigned char *out,
                               int path_len, int index_len)
{
    int total = path_len + index_len;
    memset(out, 0, (total + 3) / 4);
    node = node->pre_node; /* skip the one-base lookahead */
    for (int pos = total - 1; pos >= index_len; --pos) {
        out[pos / 4] |= (node->kmer & 3) << (6 - 2 * (pos % 4));
        node = node->pre_node;
    }
    uint64_t index = node->kmer;
    for (int pos = index_len - 1; pos >= 0; --pos) {
        out[pos / 4] |= (index & 3) << (6 - 2 * (pos % 4));
        index >>= 2;
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
    return dbgps3_load_filtered(k, p, n_thread, block_size, read_len, n_files, files, 0);
}

dbgps3_graph_t *dbgps3_load_filtered(int k, int p, int n_thread, int block_size,
                                     int read_len, int n_files, char **files,
                                     int count_bloom_shift)
{
    if (n_files < 1 || !files || k < 1 || k > 31 || p != DBGPS3_KC_BITS) return NULL;
    if (count_bloom_shift && (count_bloom_shift < p + 9 || count_bloom_shift > 40)) return NULL;
    for (int i = 0; i < n_files; ++i) if (!files[i]) return NULL;

    dbgps3_graph_t *g;
    CALLOC(g, 1);
    if (!g) return NULL;
    g->k = k;
    g->p = p;
    g->mask = (1ULL<<k*2) - 1;
    g->n_thread = n_thread;
    g->block_size = block_size;
    g->read_len = read_len > k ? read_len : k;
    g->count_bloom_shift = count_bloom_shift;

    /* Default index/payload parameters — can be overridden by caller */
    g->index_len = 16;
    g->index_byte_num = g->index_len >> 2;
    g->path_len = 140; /* default payload_len */

    if (count_bloom_shift) {
        g->_h = count_files_filtered(n_files, files, k, p, block_size, n_thread,
                                     read_len, count_bloom_shift);
        g->_source_files = (char**)calloc(n_files, sizeof(*g->_source_files));
        if (!g->_h || !g->_source_files) { dbgps3_destroy(g); return NULL; }
        for (int i = 0; i < n_files; ++i) {
            g->_source_files[i] = strdup(files[i]);
            if (!g->_source_files[i]) { dbgps3_destroy(g); return NULL; }
            ++g->_n_source_files;
        }
        return g;
    }

    fprintf(stderr, "[dbgps3_core] Counting k-mers from file 1 ......\n");
    g->_h = count_file(files[0], k, p, block_size, n_thread, read_len);
    if (!g->_h) {
        free(g);
        return NULL;
    }

    int i;
    for (i = 1; i < n_files; i++) {
        fprintf(stderr, "[dbgps3_core] Counting k-mers from file %d ......\n", i+1);
        if (!count_file2(files[i], (kc_c4x_t*)g->_h, k, p, block_size, n_thread, read_len)) {
            dbgps3_destroy(g);
            return NULL;
        }
    }

    return g;
}

int dbgps3_load_more(dbgps3_graph_t *g, int n_files, char **files)
{
    if (!g || !g->_h || n_files < 1 || !files) return -1;
    if (g->count_bloom_shift) {
        if (n_files > INT32_MAX - g->_n_source_files) return -1;
        int total = g->_n_source_files + n_files;
        char **all = (char**)calloc(total, sizeof(*all));
        if (!all) return -1;
        for (int i = 0; i < g->_n_source_files; ++i) all[i] = g->_source_files[i];
        for (int i = 0; i < n_files; ++i) all[g->_n_source_files + i] = files[i];
        dbgps3_graph_t *replacement = dbgps3_load_filtered(g->k, g->p, g->n_thread,
                g->block_size, g->read_len, total, all, g->count_bloom_shift);
        free(all);
        if (!replacement) return -1;
        if (dbgps3_set_bloom(replacement, g->bloom_bits_per_key) < 0) {
            dbgps3_destroy(replacement);
            return -1;
        }
        /* Keep query layout and codec ownership in g; retire only old data. */
        void *old_h = g->_h;
        char **old_files = g->_source_files;
        int old_n = g->_n_source_files;
        g->_h = replacement->_h;
        g->_source_files = replacement->_source_files;
        g->_n_source_files = replacement->_n_source_files;
        replacement->_h = old_h;
        replacement->_source_files = old_files;
        replacement->_n_source_files = old_n;
        dbgps3_destroy(replacement);
        return 0;
    }
    int bits_per_key = g->bloom_bits_per_key;
    c4x_set_bloom((kc_c4x_t*)g->_h, 0);
    int i;
    for (i = 0; i < n_files; i++) {
        fprintf(stderr, "[dbgps3_core] Loading additional file %d ......\n", i+1);
        if (!count_file2(files[i], (kc_c4x_t*)g->_h, g->k, g->p,
                         g->block_size, g->n_thread, g->read_len)) return -1;
    }
    return dbgps3_set_bloom(g, bits_per_key);
}

int dbgps3_set_bloom(dbgps3_graph_t *g, int bits_per_key)
{
    if (!g || c4x_set_bloom((kc_c4x_t*)g->_h, bits_per_key) < 0) return -1;
    g->bloom_bits_per_key = bits_per_key;
    return 0;
}

size_t dbgps3_bloom_bytes(const dbgps3_graph_t *g)
{
    const kc_c4x_t *h = g ? (const kc_c4x_t*)g->_h : NULL;
    return h && h->bloom ? h->bloom->n_bytes : 0;
}

void dbgps3_destroy(dbgps3_graph_t *g)
{
    if (!g) return;
    c4x_destroy((kc_c4x_t*)g->_h);
    dbgps_codec_destroy(&g->codec);
    for (int i = 0; i < g->_n_source_files; ++i) free(g->_source_files[i]);
    free(g->_source_files);
    free(g);
}

int dbgps3_set_mapping(dbgps3_graph_t *g, const char *filename,
                       char *error, size_t error_size)
{
    if (!g) return -1;
    dbgps_codec_t replacement = {0};
    if (filename && dbgps_codec_load(&replacement, filename, error, error_size) < 0) return -1;
    dbgps3_graph_t proposed = *g;
    proposed.codec = replacement;
    int index_bytes;
    if (dbgps3_validate_layout(&proposed, g->index_len, g->path_len, &index_bytes) < 0) {
        snprintf(error, error_size, "index/payload lengths are incompatible with the selected mapping");
        dbgps_codec_destroy(&replacement);
        return -1;
    }
    dbgps_codec_destroy(&g->codec);
    g->codec = replacement;
    g->index_byte_num = index_bytes;
    return 0;
}

int dbgps3_validate_layout(const dbgps3_graph_t *g, int index_len, int path_len,
                           int *index_bytes)
{
    if (!g || g->k < 1 || g->k > 31 || index_len < 0 || index_len > g->k ||
        path_len < 0 || path_len > DBGPS3_MAX_PATH_LEN) return -1;
    unsigned bases = g->codec.bases ? g->codec.bases : 1;
    unsigned bits = g->codec.bases ? g->codec.bits : 2;
    if (index_len % bases || path_len % bases ||
        (index_len / bases * bits) % 8 || (path_len / bases * bits) % 8) return -1;
    if (index_bytes) *index_bytes = index_len / bases * bits / 8;
    return 0;
}

int dbgps3_encode_index(const dbgps3_graph_t *g, uint64_t index, int index_len,
                        uint64_t *encoded)
{
    int bytes;
    if (!encoded || dbgps3_validate_layout(g, index_len, 0, &bytes) < 0 ||
        index >> (bytes * 8)) return -1;
    if (g->codec.bases)
        return dbgps_codec_encode_index(&g->codec, index, index_len, encoded);
    *encoded = dbgps3_index_encode(index, index_len);
    return 0;
}

int dbgps3_parse_index(const char *text, uint64_t *index)
{
    char *end;
    if (!text || !index || *text < '0' || *text > '9') return -1;
    errno = 0;
    *index = strtoull(text, &end, 10);
    return errno || *end ? -1 : 0;
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
    if (!g || !g->_h || !result_buf) return -1;

    const dbgps_codec_t *codec = g->codec.bases ? &g->codec : NULL;
    unsigned bases = codec ? codec->bases : 1, bits = codec ? codec->bits : 2;
    if (index_byte_num < 0 || index_byte_num > 7 || (index_byte_num * 8) % bits ||
        cut_cov < 0 || cut_cov >= DBGPS3_KC_MAX ||
        !isfinite(max_cov_ratio) || max_cov_ratio < 1.0) return -1;
    int index_len = index_byte_num * 8 / bits * bases;
    if (dbgps3_validate_layout(g, index_len, path_len, NULL) < 0) return -1;

    int k = g->k;
    uint64_t mask = g->mask;
    const kc_c4x_t *h = (const kc_c4x_t*)g->_h;
    int total_bases = path_len + index_len;
    int raw_bytes = total_bases / bases * bits / 8;

    /* Check result buffer is large enough */
    int needed = total_bases + 1;
    if (result_buf_len < needed) return -1;

    /* Allocate workspace locally — no global variables, thread-safe */
    path_node_t *all_nodes = (path_node_t *)malloc(
        DBGPS3_MAX_PATH_NUM * DBGPS3_MAX_PATH_LEN * sizeof(path_node_t));
    path_node_t **term_nodes = (path_node_t **)malloc(
        DBGPS3_MAX_PATH_NUM * sizeof(path_node_t *));
    path_node_t **next_term_nodes = (path_node_t **)malloc(
        DBGPS3_MAX_PATH_NUM * sizeof(path_node_t *));
    int ret = -1;
    if (!all_nodes || !term_nodes || !next_term_nodes) {
        free(all_nodes);
        free(term_nodes);
        free(next_term_nodes);
        return -1;
    }

    all_nodes[0].pre_node = NULL;
    all_nodes[0].kmer = start_kmer;
    all_nodes[0].cov = kmer_cov_internal(dbgps3_min_hash_key(start_kmer, k), mask, h);

    term_nodes[0] = &all_nodes[0];
    int term_node_num = 1, all_node_num = 1, p = 0;
    path_node_t *p_node;
    path_node_t next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0;
    next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    int new_term_node_num;

    while(p <= path_len && term_node_num > 0){
        int i;
        new_term_node_num = 0;
        int n_p = 0;
        while(n_p < term_node_num){
            p_node = term_nodes[n_p];
            int a_cov = p_node->cov / max_cov_ratio;
            double upper_cov = p_node->cov * max_cov_ratio;
            int b_cov = upper_cov > DBGPS3_KC_MAX ? DBGPS3_KC_MAX : (int)upper_cov;
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
                int eligible = p < skip_ratios
                    ? next_nodes[i].cov >= cut_cov + 1
                    : next_nodes[i].cov >= a_cov && next_nodes[i].cov <= b_cov;
                if (eligible) {
                    /* Stop before writing either workspace past its capacity. */
                    if (new_term_node_num >= DBGPS3_MAX_PATH_NUM ||
                        all_node_num >= DBGPS3_MAX_PATH_NUM * DBGPS3_MAX_PATH_LEN)
                        goto cleanup;
                    all_nodes[all_node_num].pre_node = p_node;
                    all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                    all_nodes[all_node_num].cov = next_nodes[i].cov;
                    next_term_nodes[new_term_node_num++] = &all_nodes[all_node_num++];
                }
            }
            n_p++;
        }
        path_node_t **swap = term_nodes;
        term_nodes = next_term_nodes;
        next_term_nodes = swap;
        term_node_num = new_term_node_num;
        ++p;
    }

    /* Full-rate codes allow every DNA word; sparse codes may reject paths.
       Uniqueness is established before touching the caller's result buffer. */
    int check_words = codec && codec->bits < 2 * codec->bases;
    if (term_node_num != 1 && !check_words) goto cleanup;
    unsigned char candidate[(DBGPS3_MAX_PATH_LEN + 31 + 3) / 4];
    unsigned char selected[sizeof(candidate)], raw[sizeof(candidate)];
    int found = 0;
    for (int i = 0; i < term_node_num; ++i) {
        if (path_len % 4 == 0 && index_len % 4 == 0)
            nodes_path(term_nodes[i], candidate, path_len, k, index_len / 4);
        else nodes_path_packed(term_nodes[i], candidate, path_len, index_len);
        if (check_words && dbgps_codec_decode(codec, candidate, total_bases, raw, raw_bytes) < 0)
            continue;
        if (found) goto cleanup;
        memcpy(selected, candidate, (total_bases + 3) / 4);
        found = 1;
    }
    if (found) {
        for (int i = 0; i < total_bases; ++i)
            result_buf[i] = dbgps3_nt4_seq_table[(selected[i / 4] >> (6 - 2 * (i % 4))) & 3];
        result_buf[total_bases] = '\0';
        ret = 0;
    }

cleanup:
    free(all_nodes);
    free(term_nodes);
    free(next_term_nodes);
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
    return index_len >= 32 ? index : index & ((1ULL << (index_len * 2)) - 1);
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
