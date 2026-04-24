#include <stdio.h>
//#include <stdlib.h>
//#include <conio.h>
//#include <curses.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop


#include "kseq.h" // FASTA/Q parser
KSEQ_INIT(gzFile, gzread)

#include "khashl.h" // hash table
#define KC_BITS 14
#define KC_MAX ((1<<KC_BITS) - 1)
#define kc_c4_eq(a, b) ((a)>>KC_BITS == (b)>>KC_BITS) // lower 10 bits for counts; higher bits for k-mer
#define kc_c4_hash(a) ((a)>>KC_BITS)

#define Max_Path_Num 1050
#define Max_Path_Len 300


KHASHL_SET_INIT(, kc_c4_t, kc_c4, uint64_t, kc_c4_hash, kc_c4_eq)

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
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
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const unsigned char nt4_seq_table[4] = { // translate 0123 to ACGT
        'A', 'C', 'G', 'T'};


static inline uint64_t hash64(uint64_t key, uint64_t mask) // invertible integer hash function
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

// The inversion of hash64(). Modified from <https://naml.us/blog/tag/invertible>
static inline uint64_t hash64i(uint64_t key, uint64_t mask) //
{
    uint64_t tmp;
    // Invert key = key + (key << 31)
    tmp = (key - (key << 31)); 	key = (key - (tmp << 31)) & mask;
    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28; 	key = key ^ tmp >> 28;
    // Invert key *= 21
    key = (key * 14933078535860113213ull) & mask;
    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14; 	tmp = key ^ tmp >> 14; 	tmp = key ^ tmp >> 14; 	key = key ^ tmp >> 14;
    // Invert key *= 265
    key = (key * 15244667743933553977ull) & mask;
    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24; 	key = key ^ tmp >> 24;
    // Invert key = (~key) + (key << 21)
    tmp = ~key; 	tmp = ~(key - (tmp << 21)); tmp = ~(key - (tmp << 21)); key = ~(key - (tmp << 21)) & mask;
    return key;
}


unsigned char* uint64_acgt(uint64_t key, unsigned char* seq, unsigned char km_len){ //decode uint_64_t kmer  to actg
    int p = km_len-1;
    while (p >= 0){
        int n = key % 4;
        key = key >> 2;
        *(seq + p) = nt4_seq_table[n];
        p = p - 1;
    }
    return seq;
}

unsigned char* uint64_int8(uint64_t key, unsigned char* seq){
    //decode uint_64_t kmer  to  char

    uint64_t mask = (1ULL<<8)-1;
    int i =0;
    for(i =0; i< 8; i++){
        *(seq+7-i) = key & mask;
        mask = mask <<8;
    }
    return seq;
}



uint64_t comp_rev2(uint64_t x, unsigned char km_len){
    uint64_t a, b=0ULL, mask = (1ULL<<km_len*2) - 1;

    a = ~x; a = a & mask;
    while((a & mask) > 0 ){
        b = b<<2;
        b = b + (a & 3ULL);
        a = a>>2;
    }
//    b = ~b;
    b = b & mask;
    return b;
}


uint64_t comp_rev(uint64_t x, unsigned char km_len){
//    if (c < 4) { // not an "N" base
//        x[0] = (x[0] << 2 | c) & mask;                  // forward strand
//        x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
//    uint64_t y=0;// mask = (1ULL << km_len * 2) - 1;
    uint64_t y=0;
    int i, c;
    for(i=0;i<km_len;i++){
        y = y << 2;
        c = x & 3ULL;
        y = y | (3-c);
        x = x >> 2;
    }
    return y;
}


uint64_t min_hash_key(uint64_t x, unsigned char km_len){
    uint64_t y = comp_rev(x, km_len);
    return (x < y) ? x:y;
}


uint64_t actgkmer_uint64(unsigned char* seq, unsigned char km_len){
//    Just for testing
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x, mask = (1ULL<<km_len*2) - 1; //shift = (km_len - 1) * 2
    for (i = 0, x = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
//            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x = 0; // if there is an "N", restart
    }

    //fprintf(stdout, "y %ld\n", x);
    return x;
}


uint64_t actgkmer_hashkey(unsigned char* seq, unsigned char km_len){
//// Function verified 20210708 Lifu Song
    int i;
    uint64_t x[2], y,  shift = (km_len - 1) * 2, mask = (1ULL<<km_len*2) - 1;
    for (i = 0, x[0] = x[1] = 0; i < km_len; ++i) {
        //fprintf(stdout, "%c ", seq[i]);
        int c = seq_nt4_table[(uint8_t)seq[i]];
        //fprintf(stdout, "%d ", c);
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
        } else x[0] = x[1] = 0; // if there is an "N", restart
    }
    y = x[0] < x[1]? x[0] : x[1];
    //fprintf(stdout, "y %ld\n", y);
    return y;
}

unsigned long uint8_actg(unsigned char bts){
    unsigned int i=0;
    unsigned long four_base_seq=0;

    for(i =0 ; i < 4; ++i){
        int m = bts & 3;
        bts = bts >>2;
        four_base_seq =four_base_seq >>8;
        four_base_seq = four_base_seq + ((unsigned long)nt4_seq_table[m]<< 24);
    }
    return four_base_seq;
}


unsigned char actg_uint8(unsigned long actg_seq){
    unsigned int i=0;
    unsigned char uint8_value=0;

    for(i =0 ; i < 4; ++i){
        unsigned int m = actg_seq & 255;
        actg_seq = actg_seq >> 8;
        uint8_value = uint8_value >> 2;
        uint8_value = uint8_value + ((unsigned char)seq_nt4_table[m]<<6);
    }
    return uint8_value;
}

void actg_intseq(unsigned char* actg_seq, unsigned char* int_seq, int8_t actg_seq_len){

    int p = 0;
    unsigned long n, n1, n2, n3, n4;
    while(p < actg_seq_len-3){
        n1 = actg_seq[p];
        n2 = actg_seq[p+1];
        n3 = actg_seq[p+2];
        n4 = actg_seq[p+3];
        n = (n1<<24) + (n2<<16) + (n3<<8) + n4;
        int_seq[p>>2] = actg_uint8(n) ;
        p = p + 4;
    }
    int_seq[p + 1] = '\0';
}


void intseq_actg(unsigned char* int_seq, unsigned char* actg_seq, int8_t int_seq_len){
//// int seq to DNA seq
//// Function not tested yet
//// one char(byte) to four bases
    int p = 0;
    unsigned long n, mask;

    unsigned char n1, n2, n3, n4;
    while(p < int_seq_len){
        mask= (1<<8) - 1;
        n = uint8_actg(int_seq[p]) ;
        n4 = (n & mask);
        mask = mask <<8;
        n3 = (n & mask)>>8;
        mask = mask <<8;
        n2 = (n & mask)>>16;
        mask = mask <<8;
        n1 = (n & mask)>>24;

        actg_seq[p<<2] = n1;
        actg_seq[(p<<2) + 1] = n2;
        actg_seq[(p<<2) + 2] = n3;
        actg_seq[(p<<2) + 3] = n4;
        p = p + 1;
    }
}

typedef struct {
    int p; // suffix length; at least 8
    kc_c4_t **h; // 1<<p hash tables
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

typedef struct {
    int n, m;
    uint64_t *a;
} buf_c4_t;


static inline void c4x_insert_buf(buf_c4_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
    int pre = y & ((1<<p) - 1);
    buf_c4_t *b = &buf[pre];
    if (b->n == b->m) {
        b->m = b->m < 8? 8 : b->m + (b->m>>1);
        REALLOC(b->a, b->m);
    }
    b->a[b->n++] = y;
}

static void count_seq_buf(buf_c4_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
    int i, l;
    uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x[0] = (x[0] << 2 | c) & mask;                  // forward strand
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
            if (++l >= k) { // we find a k-mer
                uint64_t y = x[0] < x[1]? x[0] : x[1];
                c4x_insert_buf(buf, p, hash64(y, mask));
            }
        } else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
    }
}

typedef struct { // global data structure for kt_pipeline()
    int k, block_len, n_thread, read_len;
    kseq_t *ks;
    kc_c4x_t *h;
} pldat_t;

typedef struct { // data structure for each step in kt_pipeline()
    pldat_t *p;
    int n, m, sum_len, nk;
    int *len;
    char **seq;
    buf_c4_t *buf;
} stepdat_t;

static void worker_for(void *data, long i, int tid) // callback for kt_for()
{
    stepdat_t *s = (stepdat_t*)data;
    buf_c4_t *b = &s->buf[i];
    kc_c4_t *h = s->p->h->h[i];

    int j, p = s->p->h->p;
    for (j = 0; j < b->n; ++j) {
        khint_t k;
        int absent;
        k = kc_c4_put(h, b->a[j]>>p<<KC_BITS, &absent);
        if ((kh_key(h, k)&KC_MAX) < KC_MAX) ++kh_key(h, k);
    }
}

static void *worker_pipeline(void *data, int step, void *in) // callback for kt_pipeline()
{
    pldat_t *p = (pldat_t*)data;
    if (step == 0) { // step 1: read a block of sequences
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
    } else if (step == 1) { // step 2: extract k-mers
        stepdat_t *s = (stepdat_t*)in;
        int i, n = 1<<p->h->p, m;
        CALLOC(s->buf, n);
        m = (int)(s->nk * 1.2 / n) + 1;
        for (i = 0; i < n; ++i) {
            s->buf[i].m = m;
            MALLOC(s->buf[i].a, m);
        }
        for (i = 0; i < s->n; ++i) {
            //read_len
            if(p->read_len > 0 && p->read_len < s->len[i]){
                count_seq_buf(s->buf, p->k, p->h->p, p->read_len, s->seq[i]);
            }else{
                count_seq_buf(s->buf, p->k, p->h->p, s->len[i], s->seq[i]);
            }
            free(s->seq[i]);
        }
        free(s->seq); free(s->len);
        return s;
    } else if (step == 2) { // step 3: insert k-mers to hash table
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
    if(read_len > k){
        pl.read_len = read_len;
    }else{
        pl.read_len = k;
    }

    pl.n_thread = n_thread;
    pl.h = c4x_init(p);
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}

static kc_c4x_t *count_file2(const char *fn, void *hh, int k, int p, int block_size, int n_thread, int read_len)
{
    pldat_t pl;
    gzFile fp;
    if ((fp = gzopen(fn, "r")) == 0) return 0;
    pl.ks = kseq_init(fp);
    pl.k = k;
    if(read_len > k){
        pl.read_len = read_len;
    }else{
        pl.read_len = k;
    }

    pl.n_thread = n_thread;
    pl.h = (kc_c4x_t *)hh;
    pl.block_len = block_size;
    kt_pipeline(3, worker_pipeline, &pl, 3);
    kseq_destroy(pl.ks);
    gzclose(fp);
    return pl.h;
}


typedef struct {
    uint64_t c[256];
} buf_cnt_t;

typedef struct {
    const kc_c4x_t *h;
    buf_cnt_t *cnt;
} hist_aux_t;

unsigned short kmer_cov(uint64_t kmer, uint64_t mask, kc_c4x_t *h){

    int j, x, cov=0;
    uint64_t hash_key = hash64(kmer, mask);
    j = hash_key & ((1<<KC_BITS) - 1);
    if(kh_size(h->h[j]) < 1){return 0;}

    hash_key = hash_key >> KC_BITS<< KC_BITS;
    x = kc_c4_get(h->h[j], hash_key);
    //if(x == 0){return 50;} //To avoid of segment fault
    if(kh_exist(h->h[j], x)){
        cov = kh_key(h->h[j], x) & KC_MAX;
    }
    return cov;
}

// Forward declaration to resolve circular reference
typedef struct path_node path_node;

struct path_node {
    path_node* pre_node;
    uint64_t kmer;
    unsigned int cov;
};

const uint64_t BI_MASK = 3ULL;


unsigned char four_kmers_to_byte(uint64_t k1, uint64_t k2, uint64_t k3, uint64_t k4){
    unsigned char m;
    m = (unsigned char)(((k1 & BI_MASK)<<6) + ((k2 & BI_MASK)<<4) + ((k3 & BI_MASK)<<2) + ((k4 & BI_MASK)));
    return m;
}



const uint64_t BYT_MASK = 255ULL;
unsigned char* nodes_path(path_node *term_node, unsigned char *path_seq, int path_len, int km_len, int index_byte_num){

    path_node *p_node;
    p_node = term_node;
    int path_seq_p = (path_len>>2) -1 + index_byte_num, nd_p=0;
    uint64_t kms[4];
    kms[nd_p] = p_node->kmer;

    int node_num=0;

    while(p_node->pre_node && path_seq_p >= index_byte_num){
        node_num++;
        p_node = p_node->pre_node;
        kms[nd_p] = p_node->kmer;

        if(nd_p >= 3){
            *(path_seq + path_seq_p) = four_kmers_to_byte(kms[3], kms[2], kms[1], kms[0]);
            nd_p = 0;
            path_seq_p--;
        }else{nd_p++;}
    }

    //putting index into path_seq

    while(p_node->pre_node){
        p_node = p_node->pre_node;
    }
    uint64_t index_kmer = p_node->kmer;
    int i;
    for(i=index_byte_num-1; i>=0; i--){
        *(path_seq + i) = (unsigned char) (index_kmer & BYT_MASK);
        index_kmer = index_kmer>>8;
    }


    //unsigned char DNA_seq[Max_Path_Len];

    //DNA_seq[path_len + index_byte_num*4] = '\0';
    //intseq_actg(path_seq, DNA_seq, (path_len>>2) + index_byte_num);
    //fprintf(stdout, "%s\n", DNA_seq);
    return path_seq;
}

void print_uint64_kmer(uint64_t km, int km_len){
    fprintf(stdout, "uint64: %"PRIu64, km);
   fprintf(stdout, " rev: %"PRIu64, comp_rev(km, km_len));
   fprintf(stdout, " minHash: %"PRIu64, min_hash_key(km,km_len));
    unsigned char seq[33];
    seq[32] = '\0';
    uint64_acgt(km, seq, km_len);
    fprintf(stdout, "\tkm: %s\n", seq);

}

path_node all_nodes[Max_Path_Num * Max_Path_Len];
path_node *term_nodes[Max_Path_Num+1];


int find_paths(kc_c4x_t *h, unsigned char *paths[], uint64_t kmer, uint64_t mask, int cut_cov, int path_len, int k, int index_byte_num, double max_cov_ratio, int skip_ratios){

    all_nodes[0].kmer = kmer;
    all_nodes[0].cov = kmer_cov(min_hash_key(kmer, k),mask,h);
    int a_cov = all_nodes[0].cov/max_cov_ratio, b_cov = all_nodes[0].cov * max_cov_ratio;
    if(a_cov < cut_cov + 1){a_cov = cut_cov + 1;}

    term_nodes[0] = &all_nodes[0];

    int term_node_num = 1, all_node_num = 1, p = 0;
    int seq_len = (path_len>>2) + index_byte_num;

    path_node *p_node;
    path_node next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0; next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    int new_term_node_num;

    while(p <= path_len && term_node_num > 0 && term_node_num <= Max_Path_Num){
        int i;

        new_term_node_num = 0;

        int n_p = 0;
        while(n_p < term_node_num && new_term_node_num <=Max_Path_Num){

            p_node = term_nodes[n_p];

            a_cov = p_node->cov/max_cov_ratio, b_cov = p_node->cov * max_cov_ratio;
            if(a_cov < cut_cov + 1){a_cov = cut_cov + 1;}

            next_nodes[0].kmer = ((*p_node).kmer<<2) & mask;
            next_nodes[1].kmer = (((*p_node).kmer<<2) + 1) & mask;
            next_nodes[2].kmer = (((*p_node).kmer<<2) + 2) & mask;
            next_nodes[3].kmer = (((*p_node).kmer<<2) + 3) & mask;

            next_nodes[0].cov = kmer_cov(min_hash_key(next_nodes[0].kmer, k),mask,h);
            next_nodes[1].cov = kmer_cov(min_hash_key(next_nodes[1].kmer, k),mask,h);
            next_nodes[2].cov = kmer_cov(min_hash_key(next_nodes[2].kmer, k),mask,h);
            next_nodes[3].cov = kmer_cov(min_hash_key(next_nodes[3].kmer, k),mask,h);

            // Reuse the same i from the outer loop
            for(i=0;i<4;i++){
                // Check if we're still in the skipping phase (skip ratio checks for first skip_ratios nodes)
                if (p < skip_ratios) {
                    // During skipping phase, only check minimum coverage
                    if(next_nodes[i].cov >= cut_cov + 1) {
                        all_nodes[all_node_num].pre_node = p_node;
                        all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                        all_nodes[all_node_num].cov = next_nodes[i].cov;
                        term_nodes[new_term_node_num] = &all_nodes[all_node_num];
                        new_term_node_num++;
                        all_node_num++;
                    }
                } else {
                    // Normal phase: check ratio bounds as before
                    if(next_nodes[i].cov >=a_cov && next_nodes[i].cov <= b_cov) {

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

    if(term_node_num == 1){
        CALLOC(paths[0], seq_len);
        nodes_path(term_nodes[0], paths[0], path_len, k, index_byte_num);
        return 0;
    }else{return -1;}
}


uint64_t primer_uint64( char *pr, int pr_len, int index_len, int km_len){
    unsigned char *seq;
    uint64_t pr_int64;
    int i = pr_len + index_len - km_len, j=0;

    MALLOC(seq, km_len);

    while( i < pr_len ){
        *(seq + j) = *(pr + i);
        i++;
        j++;
    }

    while(j < km_len){
        *(seq + j) = 'A';
        j++;
    }

    pr_int64 = actgkmer_uint64(seq, km_len);
    free(seq);
    return pr_int64;
}

void print_kmer_seq(unsigned char *seq, int k){ //Debugging function
    int i=0;
    for(i=0;i<k;i++){
            printf("%c", *(seq + i));
    }
}

void print_cov(unsigned char* seq, int k, kc_c4x_t *h){ //Debugging function
    uint64_t mask = (1ULL<<k*2) - 1;
    int test_cov=0;
    test_cov = kmer_cov( actgkmer_hashkey(seq, k), mask, h);
    print_kmer_seq(seq, k);
    printf("\t");
    printf("cov: %d\n", test_cov);
}


int main(int argc, char *argv[])
{
    kc_c4x_t *h;
    int i, c=1, k = 31, p = KC_BITS, block_size = 10000000, n_thread = 3, min_cov = 1, index_len = 16, payload_len=140, path_len;
    double max_cov_ratio = 2.55; // Default value for Max_Cov_Ratio
    double max_r = 0 ; // Maximum r value for ratio range iteration, default 0 (disabled)
    int mem_test = 0;
    int read_len = 200;
    int skip_ratios = 0; // Number of initial ratio values to skip, default 0
    char *output_prefix = NULL; // Output prefix for log file

    unsigned int min_index = 101010102, max_index=101295684, index_byte_num, max_cov=0;
    uint64_t mask, pr_uint64;

    char p1[] = "CCTGCAGAGTAGCATGTC"; //p2[18] = 'CTGACACTGATGCATCCG';
    int pr_len = 18;


    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:c:l:a:b:i:C:r:R:m:s:o:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'C') max_cov = atoi(o.arg);

        else if (c == 'a') min_index = atoi(o.arg);
        else if (c == 'b') max_index = atoi(o.arg);
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 'r') max_cov_ratio = atof(o.arg); // Changed this to use max_cov_ratio
        else if (c == 'R') max_r = atof(o.arg); // Maximum r value for ratio range iteration
        else if (c == 'm') mem_test = 1;
        else if (c == 's') skip_ratios = atoi(o.arg); // Number of initial ratio values to skip
        else if (c == 'o') output_prefix = o.arg; // Output prefix for log file
    }

    if (argc - o.ind < 1) {
        fprintf(stderr, "\n************************************************************************\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "**  DBGPS v2 - De Bruijn Graph based strand assembller for DNA data storage  **\n");
        fprintf(stderr, "**     Version 2.0 (2025-05-14)  Author: Lifu Song lifu.song@outlook.com      **\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "************************************************************************\n\n");
        fprintf(stderr, "Usage: DBGPS [options] <input file> \n");
        fprintf(stderr, "                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
        fprintf(stderr, "  -i INT     length of index [%d] bp\n", index_len);
        fprintf(stderr, "  -l INT     data encoding length [%d] bp\n", payload_len);
        //fprintf(stderr, "  -p INT     prefix length [%d]\n", p);
        fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
        fprintf(stderr, "  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [%d]\n", min_cov);
        fprintf(stderr, "  -C INT     Switch on k-mer coverage testing mode.  [%d]\n", max_cov);

        fprintf(stderr, "  -a INT     Initial index [%d]\n", min_index);
        fprintf(stderr, "  -b INT     End index [%d]\n", max_index);

        fprintf(stderr, "  -r FLOAT   Maximum coverage ratio [%.1f]\n", max_cov_ratio);
        fprintf(stderr, "  -R FLOAT   Upper bound for ratio range iteration (0 = disabled) [%.1f]\n", max_r);
        fprintf(stderr, "  -s INT     Number of initial positions to skip ratio checks for [%d]\n", skip_ratios);
        fprintf(stderr, "  -o STR     Output prefix for log file [none]\n");

        fprintf(stderr, "\n");
        return 1;
    }
    if (p < KC_BITS) {
        fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
        return 1;
    }


    if(max_cov < min_cov){max_cov = min_cov;}
    index_byte_num = index_len>>2;
    mask = (1ULL<<k*2) - 1;
    pr_uint64 = primer_uint64(p1, pr_len, index_len, k);

    path_len = payload_len;

    clock_t start,end;
    double kc_dur, dec_dur;
    
    // Files for output if output_prefix is provided
    FILE *log_file = NULL;
    FILE *strands_file = NULL;
    if (output_prefix != NULL) {
        // Open log file
        char log_filename[1024];
        sprintf(log_filename, "%s.log", output_prefix);
        log_file = fopen(log_filename, "w");
        if (log_file == NULL) {
            fprintf(stderr, "Warning: Could not open log file %s for writing\n", log_filename);
        } else {
            fprintf(stderr, "Writing parameters to log file: %s\n", log_filename);
        }
        
        // Open strands file
        char strands_filename[1024];
        sprintf(strands_filename, "%s.strands", output_prefix);
        strands_file = fopen(strands_filename, "w");
        if (strands_file == NULL) {
            fprintf(stderr, "Warning: Could not open strands file %s for writing\n", strands_filename);
        } else {
            fprintf(stderr, "Will write assembled strands to: %s\n", strands_filename);
        }
    }
    
    // Function to print to both stdout and log file
    #define PRINT_BOTH(format, ...) do { \
        fprintf(stdout, format, ##__VA_ARGS__); \
        if (log_file != NULL) fprintf(log_file, format, ##__VA_ARGS__); \
    } while(0)
    
    PRINT_BOTH("========== DBGPS Parameters ==========\n");
    PRINT_BOTH("k-mer size (-k): %d\n", k);
    PRINT_BOTH("Prefix length (-p): %d\n", p);
    PRINT_BOTH("Number of threads (-t): %d\n", n_thread);
    PRINT_BOTH("Minimum coverage cutoff (-c): %d\n", min_cov);
    PRINT_BOTH("Maximum coverage for testing (-C): %d\n", max_cov);
    PRINT_BOTH("Index length (-i): %d bp\n", index_len);
    PRINT_BOTH("Data encoding length (-l): %d bp\n", payload_len);
    PRINT_BOTH("Read length for k-mer counting: %d\n", read_len);
    PRINT_BOTH("Skip ratio checks for first (-s): %d positions\n", skip_ratios);
    PRINT_BOTH("Initial index (-a): %d\n", min_index);
    PRINT_BOTH("End index (-b): %d\n", max_index);
    PRINT_BOTH("Maximum coverage ratio (-r): %.2f\n", max_cov_ratio);
    PRINT_BOTH("Ratio range upper bound (-R): %.2f\n", max_r);
    PRINT_BOTH("Memory test mode: %s\n", mem_test ? "enabled" : "disabled");
    PRINT_BOTH("Path length (payload): %d\n", path_len);
    PRINT_BOTH("Index byte number: %d\n", index_byte_num);
    PRINT_BOTH("Output prefix: %s\n", output_prefix ? output_prefix : "none");
    PRINT_BOTH("======================================\n");
    
    fprintf(stderr, "Counting k-mers ......\n");


    start = clock();
    fprintf(stderr, "Counting NGS file 1 ......\n");
    h = count_file(argv[o.ind], k, p, block_size, n_thread, read_len);

    int f_num = argc - o.ind, c_f_n = 2;
    while(c_f_n <= f_num ) {
        fprintf(stderr, "Counting NGS file %d ......\n", c_f_n);
        h = count_file2(argv[o.ind + c_f_n - 1], h, k, p, block_size, n_thread, read_len);  
        c_f_n = c_f_n + 1;
    }


    end = clock();
    kc_dur = (double)(end - start);

    PRINT_BOTH("\nK-mer counting time: %f seconds\n",(kc_dur/CLOCKS_PER_SEC));
    if(mem_test){
            fprintf(stderr, "Memory test mode activated! Press ENTER to continue...\n");
            getchar();}
    fprintf(stderr, "Decoding strands from index %d to ",min_index);
    fprintf(stderr, "%d ......\n",max_index);

    unsigned int id;
    unsigned int cov_cut;
    int c_index;
    int dec_strand_num=0;
    
    // Check if ratio range iteration is enabled
    if (max_r <= 0 || max_r <= max_cov_ratio) {
        // Original behavior - use only the specified max_cov_ratio
        for(cov_cut=min_cov; cov_cut<=max_cov; cov_cut++) {
            start = clock();
            dec_strand_num = 0;

            for(id=min_index; id<=max_index; id++) {
                unsigned char *path[Max_Path_Num];
                c_index = find_paths(h, path, pr_uint64 | (uint64_t)id, mask, cov_cut, path_len, k, index_byte_num, max_cov_ratio, skip_ratios);
                if(c_index > -1) {
                    dec_strand_num++;
                    unsigned char DNA_seq[Max_Path_Len];

                    DNA_seq[path_len + index_byte_num*4] = '\0';
                    intseq_actg(path[c_index], DNA_seq, (path_len>>2) + index_byte_num);
                    
                    if (output_prefix != NULL && strands_file != NULL) {
                        // Output to strands file if -o is specified
                        fprintf(strands_file, "%d\t%s\n", id, DNA_seq);
                    } else {
                        // Output to stdout if no output file is specified
                        printf("%d\t%s\n", id, DNA_seq);
                    }
                }
            }

            end = clock();
            dec_dur = (double)(end - start);
            // Update to use PRINT_BOTH
            PRINT_BOTH("\nCoverage cutoff: %d, Ratio: %.2f", cov_cut, max_cov_ratio);
            PRINT_BOTH("\nStrands decoded: %d ", dec_strand_num);
            PRINT_BOTH("\nStrand decoding time: %f seconds\n",(dec_dur/CLOCKS_PER_SEC));
        }
    } else {
        // New behavior - iterate through different r values
        for(double r = max_cov_ratio; r <= max_r; r += 0.05) {
            for(cov_cut=min_cov; cov_cut<=max_cov; cov_cut++) {
                start = clock();
                dec_strand_num = 0;

                for(id=min_index; id<=max_index; id++) {
                    unsigned char *path[Max_Path_Num];
                    c_index = find_paths(h, path, pr_uint64 | (uint64_t)id, mask, cov_cut, path_len, k, index_byte_num, r, skip_ratios);
                    if(c_index > -1) {
                        dec_strand_num++;
                        unsigned char DNA_seq[Max_Path_Len];

                        DNA_seq[path_len + index_byte_num*4] = '\0';
                        intseq_actg(path[c_index], DNA_seq, (path_len>>2) + index_byte_num);
                        
                        if (output_prefix != NULL && strands_file != NULL) {
                            // Output to strands file if -o is specified
                            fprintf(strands_file, "%d\t%s\n", id, DNA_seq);
                        } else {
                            // Output to stdout if no output file is specified
                            printf("%d\t%s\n", id, DNA_seq);
                        }
                    }
                }

                end = clock();
                dec_dur = (double)(end - start);
                // Update to use PRINT_BOTH
                PRINT_BOTH("\nCoverage cutoff: %d, Ratio: %.2f", cov_cut, r);
                PRINT_BOTH("\nStrands decoded: %d ", dec_strand_num);
                PRINT_BOTH("\nStrand decoding time: %f seconds\n",(dec_dur/CLOCKS_PER_SEC));
            }
        }
    }

    for (i = 0; i < 1<<p; ++i) kc_c4_destroy(h->h[i]);
    free(h->h); free(h);
    
    // Close output files if they were opened
    if (log_file != NULL) {
        fclose(log_file);
    }
    
    if (strands_file != NULL) {
        fclose(strands_file);
    }
    
    #undef PRINT_BOTH
    return 0;
}




/*

The CRC function was modified from a python implementation:


"""Pure python library for calculating CRC16"""

##############################################################################
#
#    Copyright (C) Gennady Trafimenkov, 2011
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################


# table for calculating CRC
# this particular table was generated using pycrc v0.7.6, http://www.tty1.net/pycrc/
# using the configuration:
#  *    Width        = 16
#  *    Poly         = 0x1021
#  *    XorIn        = 0x0000
#  *    ReflectIn    = False
#  *    XorOut       = 0x0000
#  *    ReflectOut   = False
#  *    Algorithm    = table-driven
# by following command:
#   python pycrc.py --model xmodem --algorithm table-driven --generate c

*/
