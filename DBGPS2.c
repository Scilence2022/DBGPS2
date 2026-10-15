#include <stdio.h>
//#include <stdlib.h>
//#include <conio.h>
//#include <curses.h>
#include <stdint.h>
#include <zlib.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include "ketopt.h" // command-line argument parser
#include "kthread.h" // multi-threading models: pipeline and multi-threaded for loop
#include "dbgps_codec.h"
#include "dbgps_bloom.h"


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

#define DBGPS_COUNT_BITS KC_BITS
#define DBGPS_COUNT_MAX KC_MAX
#include "dbgps_bloom_graph.h"

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
    int read_error;
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
    for (int j = 0; j < b->n; ++j)
        c4x_count_key(s->p->h, (int)i, b->a[j]);
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
        if (ret < -1) p->read_error = 1;
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
    pldat_t pl = {0};
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
    if(read_len > k){
        pl.read_len = read_len;
    }else{
        pl.read_len = k;
    }

    pl.n_thread = n_thread;
    pl.h = (kc_c4x_t *)hh;
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


typedef struct {
    uint64_t c[256];
} buf_cnt_t;

typedef struct {
    const kc_c4x_t *h;
    buf_cnt_t *cnt;
} hist_aux_t;




static const unsigned short crc16tab[] = {
        0x0000, 0x1021, 0x2042, 0x3063, 0x4084, 0x50a5, 0x60c6, 0x70e7,
        0x8108, 0x9129, 0xa14a, 0xb16b, 0xc18c, 0xd1ad, 0xe1ce, 0xf1ef,
        0x1231, 0x0210, 0x3273, 0x2252, 0x52b5, 0x4294, 0x72f7, 0x62d6,
        0x9339, 0x8318, 0xb37b, 0xa35a, 0xd3bd, 0xc39c, 0xf3ff, 0xe3de,
        0x2462, 0x3443, 0x0420, 0x1401, 0x64e6, 0x74c7, 0x44a4, 0x5485,
        0xa56a, 0xb54b, 0x8528, 0x9509, 0xe5ee, 0xf5cf, 0xc5ac, 0xd58d,
        0x3653, 0x2672, 0x1611, 0x0630, 0x76d7, 0x66f6, 0x5695, 0x46b4,
        0xb75b, 0xa77a, 0x9719, 0x8738, 0xf7df, 0xe7fe, 0xd79d, 0xc7bc,
        0x48c4, 0x58e5, 0x6886, 0x78a7, 0x0840, 0x1861, 0x2802, 0x3823,
        0xc9cc, 0xd9ed, 0xe98e, 0xf9af, 0x8948, 0x9969, 0xa90a, 0xb92b,
        0x5af5, 0x4ad4, 0x7ab7, 0x6a96, 0x1a71, 0x0a50, 0x3a33, 0x2a12,
        0xdbfd, 0xcbdc, 0xfbbf, 0xeb9e, 0x9b79, 0x8b58, 0xbb3b, 0xab1a,
        0x6ca6, 0x7c87, 0x4ce4, 0x5cc5, 0x2c22, 0x3c03, 0x0c60, 0x1c41,
        0xedae, 0xfd8f, 0xcdec, 0xddcd, 0xad2a, 0xbd0b, 0x8d68, 0x9d49,
        0x7e97, 0x6eb6, 0x5ed5, 0x4ef4, 0x3e13, 0x2e32, 0x1e51, 0x0e70,
        0xff9f, 0xefbe, 0xdfdd, 0xcffc, 0xbf1b, 0xaf3a, 0x9f59, 0x8f78,
        0x9188, 0x81a9, 0xb1ca, 0xa1eb, 0xd10c, 0xc12d, 0xf14e, 0xe16f,
        0x1080, 0x00a1, 0x30c2, 0x20e3, 0x5004, 0x4025, 0x7046, 0x6067,
        0x83b9, 0x9398, 0xa3fb, 0xb3da, 0xc33d, 0xd31c, 0xe37f, 0xf35e,
        0x02b1, 0x1290, 0x22f3, 0x32d2, 0x4235, 0x5214, 0x6277, 0x7256,
        0xb5ea, 0xa5cb, 0x95a8, 0x8589, 0xf56e, 0xe54f, 0xd52c, 0xc50d,
        0x34e2, 0x24c3, 0x14a0, 0x0481, 0x7466, 0x6447, 0x5424, 0x4405,
        0xa7db, 0xb7fa, 0x8799, 0x97b8, 0xe75f, 0xf77e, 0xc71d, 0xd73c,
        0x26d3, 0x36f2, 0x0691, 0x16b0, 0x6657, 0x7676, 0x4615, 0x5634,
        0xd94c, 0xc96d, 0xf90e, 0xe92f, 0x99c8, 0x89e9, 0xb98a, 0xa9ab,
        0x5844, 0x4865, 0x7806, 0x6827, 0x18c0, 0x08e1, 0x3882, 0x28a3,
        0xcb7d, 0xdb5c, 0xeb3f, 0xfb1e, 0x8bf9, 0x9bd8, 0xabbb, 0xbb9a,
        0x4a75, 0x5a54, 0x6a37, 0x7a16, 0x0af1, 0x1ad0, 0x2ab3, 0x3a92,
        0xfd2e, 0xed0f, 0xdd6c, 0xcd4d, 0xbdaa, 0xad8b, 0x9de8, 0x8dc9,
        0x7c26, 0x6c07, 0x5c64, 0x4c45, 0x3ca2, 0x2c83, 0x1ce0, 0x0cc1,
        0xef1f, 0xff3e, 0xcf5d, 0xdf7c, 0xaf9b, 0xbfba, 0x8fd9, 0x9ff8,
        0x6e17, 0x7e36, 0x4e55, 0x5e74, 0x2e93, 0x3eb2, 0x0ed1, 0x1ef0,
};


unsigned short  crc16xmodem(unsigned char *data, int len){
    unsigned short crc = 0;
    int i;

    for (i = 0; i < len; ++i) {
        crc = ((crc<<8)&0xff00) ^ crc16tab[((crc>>8)&0xff)^data[i]];
    }
    return crc;
}

int crc16_check(unsigned char *data, int data_len){
     unsigned short crc = 0;
     crc = *(data + data_len -2) * 256 +  *(data + data_len -1);
     //crc = (unsigned short)*(data + data_len -2);

     unsigned short crc_cl = crc16xmodem(data, data_len-2);

     if(crc == crc_cl){
        return 1;
     }else{
        return 0;
     }
}


unsigned short kmer_cov(uint64_t kmer, uint64_t mask, kc_c4x_t *h){

    int j, x, cov=0;
    uint64_t hash_key = hash64(kmer, mask);
    if (h->bloom && !dbgps_bloom_contains(h->bloom, hash_key)) return 0;
    j = hash_key & ((1<<KC_BITS) - 1);
    if(kh_size(h->h[j]) < 1){return 0;}

    hash_key = hash_key >> KC_BITS<< KC_BITS;
    x = kc_c4_get(h->h[j], hash_key);
    //if(x == 0){return 50;} //To avoid of segment fault
    if(x != kh_end(h->h[j]) && kh_exist(h->h[j], x)){
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

    while(p_node->pre_node && path_seq_p >= index_byte_num){
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
path_node *term_nodes[Max_Path_Num];
path_node *next_term_nodes[Max_Path_Num];


/* Reconstruct arbitrary block boundaries, retaining canonical packed DNA for
   output. The terminal node is the historical one-base primer lookahead. */
static void nodes_path_packed(path_node *node, unsigned char *out, int path_len, int index_len)
{
    int total = path_len + index_len;
    memset(out, 0, (total + 3) / 4);
    node = node->pre_node;
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

static int codec_path_sizes(const dbgps_codec_t *codec, int path_len, int index_len,
                            int k, int *packed_bytes, int *raw_bytes)
{
    unsigned bases = codec ? codec->bases : 1, bits = codec ? codec->bits : 2;
    if (k < 1 || k > 31 || path_len < 0 || index_len < 0 || index_len > k ||
        path_len >= Max_Path_Len || index_len + path_len >= Max_Path_Len ||
        !bases || index_len % bases || path_len % bases) return -1;
    unsigned padding = (bits - 16 % bits) % bits;
    unsigned total_bits = (index_len + path_len) / bases * bits;
    if (index_len / bases * bits % 8 || total_bits < 16 + padding ||
        (total_bits - padding) % 8) return -1;
    *packed_bytes = (index_len + path_len + 3) / 4;
    *raw_bytes = (total_bits - padding) / 8;
    return 0;
}

static int find_paths_mapped(kc_c4x_t *h, unsigned char *paths[], uint64_t kmer, uint64_t mask, int cut_cov, int path_len, int k, int index_len, double max_cov_ratio, int skip_ratios, const dbgps_codec_t *codec){

    int seq_len, raw_len;
    if (codec_path_sizes(codec, path_len, index_len, k, &seq_len, &raw_len) < 0 ||
        cut_cov < 0 || cut_cov >= KC_MAX ||
        !isfinite(max_cov_ratio) || max_cov_ratio < 1.0) return -1;

    all_nodes[0].pre_node = NULL;
    all_nodes[0].kmer = kmer;
    all_nodes[0].cov = kmer_cov(min_hash_key(kmer, k),mask,h);

    term_nodes[0] = &all_nodes[0];

    int term_node_num = 1, all_node_num = 1, p = 0;

    path_node *p_node;
    path_node next_nodes[4];
    next_nodes[0].cov = 0; next_nodes[1].cov = 0; next_nodes[2].cov = 0; next_nodes[3].cov = 0;

    int new_term_node_num;

    while(p <= path_len && term_node_num > 0){
        int i;

        new_term_node_num = 0;

        int n_p = 0;
        while(n_p < term_node_num){

            p_node = term_nodes[n_p];

            int a_cov = p_node->cov/max_cov_ratio;
            double upper_cov = p_node->cov * max_cov_ratio;
            int b_cov = upper_cov > KC_MAX ? KC_MAX : (int)upper_cov;
            if(a_cov < cut_cov + 1){a_cov = cut_cov + 1;}

            next_nodes[0].kmer = ((*p_node).kmer<<2) & mask;
            next_nodes[1].kmer = (((*p_node).kmer<<2) + 1) & mask;
            next_nodes[2].kmer = (((*p_node).kmer<<2) + 2) & mask;
            next_nodes[3].kmer = (((*p_node).kmer<<2) + 3) & mask;

            next_nodes[0].cov = kmer_cov(min_hash_key(next_nodes[0].kmer, k),mask,h);
            next_nodes[1].cov = kmer_cov(min_hash_key(next_nodes[1].kmer, k),mask,h);
            next_nodes[2].cov = kmer_cov(min_hash_key(next_nodes[2].kmer, k),mask,h);
            next_nodes[3].cov = kmer_cov(min_hash_key(next_nodes[3].kmer, k),mask,h);

            for(i=0;i<4;i++){
                int eligible = p < skip_ratios
                    ? next_nodes[i].cov >= cut_cov + 1
                    : next_nodes[i].cov >= a_cov && next_nodes[i].cov <= b_cov;
                if (eligible) {
                    /* An incomplete search cannot establish a unique CRC-valid path. */
                    if (new_term_node_num >= Max_Path_Num ||
                        all_node_num >= Max_Path_Num * Max_Path_Len) return -1;
                    all_nodes[all_node_num].pre_node = p_node;
                    all_nodes[all_node_num].kmer = next_nodes[i].kmer;
                    all_nodes[all_node_num].cov = next_nodes[i].cov;
                    next_term_nodes[new_term_node_num++] = &all_nodes[all_node_num++];
                }
            }



            n_p++;
        }
        /* Keep this generation intact until every parent has been visited. */
        memcpy(term_nodes, next_term_nodes, new_term_node_num * sizeof(*term_nodes));
        term_node_num = new_term_node_num;
        ++p;
    }

    if (term_node_num == 0) return -1;

    /* Scratch space is bounded by Max_Path_Len; avoid two heap allocations
       for every index, cutoff, and ratio that reaches candidate validation. */
    unsigned char candidate[(Max_Path_Len + 3) / 4];
    unsigned char path_tmp[(Max_Path_Len + 3) / 4];

    int c_index = -1;
    for (int i = 0; i < term_node_num; i++) {
        unsigned char *crc_data = candidate;
        if (!codec || codec->byte_fast) {
            nodes_path(term_nodes[i], candidate, path_len, k, index_len / 4);
            if (codec) {
                for (int j = 0; j < seq_len; j++) path_tmp[j] = codec->decode_byte[candidate[j]];
                crc_data = path_tmp;
            }
        } else {
            nodes_path_packed(term_nodes[i], candidate, path_len, index_len);
            if (dbgps_codec_decode(codec, candidate, index_len + path_len, path_tmp, raw_len) < 0)
                continue;
            crc_data = path_tmp;
        }
        if (crc16_check(crc_data, raw_len)) {
            if (c_index >= 0) {
                free(paths[c_index]);
                paths[c_index] = NULL;
                c_index = -1;
                break;
            }
            paths[i] = (unsigned char *)malloc(seq_len);
            if (!paths[i]) break;
            memcpy(paths[i], candidate, seq_len);
            c_index = i;
        }
    }
    /* Only the unique result is retained; the caller owns that allocation. */
    return c_index;
}

int find_paths(kc_c4x_t *h, unsigned char *paths[], uint64_t kmer, uint64_t mask, int cut_cov, int path_len, int k, int index_byte_num, double max_cov_ratio, int skip_ratios)
{
    if (index_byte_num < 0 || index_byte_num > k / 4) return -1;
    return find_paths_mapped(h, paths, kmer, mask, cut_cov, path_len, k,
                             index_byte_num * 4, max_cov_ratio, skip_ratios, NULL);
}


uint64_t index_encode(uint64_t index, unsigned int index_len){
    /* Index bits already are canonical packed DNA: A=00,C=01,G=10,T=11. */
    return index_len >= 32 ? index : index & ((1ULL << (index_len * 2)) - 1);
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

static int parse_index(const char *text, uint64_t *value)
{
    char *end;
    if (*text < '0' || *text > '9') return -1;
    errno = 0;
    *value = strtoull(text, &end, 10);
    return errno || *end ? -1 : 0;
}

static void packed_dna_string(const unsigned char *packed, unsigned char *dna, int bases)
{
    for (int i = 0; i < bases; ++i)
        dna[i] = nt4_seq_table[(packed[i / 4] >> (6 - 2 * (i % 4))) & 3];
    dna[bases] = '\0';
}


int main(int argc, char *argv[])
{
    kc_c4x_t *h = NULL;
    int exit_status = 0;
    int c=1, k = 31, p = KC_BITS, block_size = 10000000, n_thread = 3, min_cov = 1, index_len = 16, payload_len=140, crc_len=8, path_len;
    double max_cov_ratio = 2.55; // Default value for Max_Cov_Ratio
    double max_r = 0 ; // Maximum r value for ratio range iteration, default 0 (disabled)
    int mem_test = 0;
    int bloom_bits_per_key = 0;
    int count_bloom_shift = 0;
    int read_len = 200;
    int skip_ratios = 0; // Number of initial ratio values to skip, default 0
    char *output_prefix = NULL; // Output prefix for log file
    const char *mapping_file = NULL;
    dbgps_codec_t mapping = {0};
    const dbgps_codec_t *codec = NULL; /* NULL uses canonical packed DNA directly. */

    uint64_t min_index = 101010102, max_index = 101295684;
    unsigned int max_cov = 0;
    uint64_t mask, pr_uint64;

    char p1[] = "CCTGCAGAGTAGCATGTC"; //p2[18] = 'CTGACACTGATGCATCCG';
    int pr_len = 18;


    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:c:l:a:b:i:C:r:R:m:s:o:M:B:F:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'C') max_cov = atoi(o.arg);

        else if (c == 'a' || c == 'b') {
            if (parse_index(o.arg, c == 'a' ? &min_index : &max_index) < 0) {
                fprintf(stderr, "ERROR: -%c requires an unsigned decimal index\n", c);
                return 1;
            }
        }
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 'r') max_cov_ratio = atof(o.arg); // Changed this to use max_cov_ratio
        else if (c == 'R') max_r = atof(o.arg); // Maximum r value for ratio range iteration
        else if (c == 'm') mem_test = 1;
        else if (c == 's') skip_ratios = atoi(o.arg); // Number of initial ratio values to skip
        else if (c == 'o') output_prefix = o.arg; // Output prefix for log file
        else if (c == 'M') mapping_file = o.arg;
        else if (c == 'F') {
            if (dbgps_bloom_parse_shift(o.arg, &count_bloom_shift) < 0) {
                fprintf(stderr, "ERROR: -F requires 0 (disabled) or log2 bitmap bits from 23 to 40\n");
                return 1;
            }
        }
        else if (c == 'B') {
            if (dbgps_bloom_parse_bits(o.arg, &bloom_bits_per_key) < 0) {
                fprintf(stderr, "ERROR: -B requires an integer from 0 to 64 (0 disables Bloom filtering)\n");
                return 1;
            }
        }
        else {
            fprintf(stderr, "ERROR: unknown option or missing option value\n");
            return 1;
        }
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

        fprintf(stderr, "  -a INT     Initial index [%"PRIu64"]\n", min_index);
        fprintf(stderr, "  -b INT     End index [%"PRIu64"]\n", max_index);

        fprintf(stderr, "  -r FLOAT   Maximum coverage ratio [%.1f]\n", max_cov_ratio);
        fprintf(stderr, "  -R FLOAT   Upper bound for ratio range iteration (0 = disabled) [%.1f]\n", max_r);
        fprintf(stderr, "  -s INT     Number of initial positions to skip ratio checks for [%d]\n", skip_ratios);
        fprintf(stderr, "  -o STR     Output prefix for log file [none]\n");
        fprintf(stderr, "  -M FILE    Custom DNA/bit mapping [none; built-in A=00,C=01,G=10,T=11]\n");
        fprintf(stderr, "  -B INT     Bloom query prefilter bits/key (1..64; 0 disables) [0]\n");
        fprintf(stderr, "  -F INT     Two-pass singleton filter log2 bits (23..40; 0 disables) [0]\n");

        fprintf(stderr, "\n");
        return 1;
    }
    if (p < KC_BITS) {
        fprintf(stderr, "ERROR: -p should be at least %d\n", KC_BITS);
        return 1;
    }

    if (mapping_file) {
        char error[512];
        if (dbgps_codec_load(&mapping, mapping_file, error, sizeof(error)) < 0) {
            fprintf(stderr, "ERROR: %s\n", error);
            return 1;
        }
        codec = &mapping;
        crc_len = ((16 + codec->bits - 1) / codec->bits) * codec->bases;
    }

    if(max_cov < min_cov){max_cov = min_cov;}
    int packed_bytes, raw_bytes;
    unsigned bases = codec ? codec->bases : 1, bits = codec ? codec->bits : 2;
    if (payload_len < 0 || payload_len >= Max_Path_Len || payload_len % bases ||
        (payload_len / bases * bits) % 8 ||
        codec_path_sizes(codec, payload_len + crc_len, index_len, k, &packed_bytes, &raw_bytes) < 0 ||
        k - index_len > pr_len || min_index > max_index ||
        max_index >> (index_len / bases * bits)) {
        fprintf(stderr, "ERROR: incompatible strand layout or index range; index and payload must\n"
                        "       contain whole mapping blocks and whole bytes, index must fit k,\n"
                        "       k-index must fit the primer, and total strand length must be < %d bp\n", Max_Path_Len);
        dbgps_codec_destroy(&mapping);
        return 1;
    }
    mask = (1ULL<<k*2) - 1;
    pr_uint64 = primer_uint64(p1, pr_len, index_len, k);

    path_len = payload_len + crc_len;
    if (codec) read_len = 2 * pr_len + index_len + path_len;

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
    PRINT_BOTH("Bloom counting filter (-F): %d (0 disables; singletons removed when enabled)\n", count_bloom_shift);
    PRINT_BOTH("Minimum coverage cutoff (-c): %d\n", min_cov);
    PRINT_BOTH("Maximum coverage for testing (-C): %d\n", max_cov);
    PRINT_BOTH("Index length (-i): %d bp\n", index_len);
    PRINT_BOTH("Data encoding length (-l): %d bp\n", payload_len);
    PRINT_BOTH("Read length for k-mer counting: %d\n", read_len);
    PRINT_BOTH("Skip ratio checks for first (-s): %d positions\n", skip_ratios);
    PRINT_BOTH("Initial index (-a): %"PRIu64"\n", min_index);
    PRINT_BOTH("End index (-b): %"PRIu64"\n", max_index);
    PRINT_BOTH("Maximum coverage ratio (-r): %.2f\n", max_cov_ratio);
    PRINT_BOTH("Ratio range upper bound (-R): %.2f\n", max_r);
    PRINT_BOTH("Memory test mode: %s\n", mem_test ? "enabled" : "disabled");
    PRINT_BOTH("Path length (payload + CRC): %d\n", path_len);
    PRINT_BOTH("Index byte number: %u\n", index_len / bases * bits / 8);
    if (codec) {
        PRINT_BOTH("Mapping: %s (%u bases <-> %u bits; %s)\n", mapping_file,
                   bases, bits, codec->byte_fast ? "byte lookup" : "block lookup");
    } else {
        PRINT_BOTH("Mapping: none (built-in A=00,C=01,G=10,T=11; direct 2-bit packing)\n");
    }
    PRINT_BOTH("CRC encoding length: %d bp; zero padding: %u bits\n", crc_len, (bits - 16 % bits) % bits);
    PRINT_BOTH("Output prefix: %s\n", output_prefix ? output_prefix : "none");
    PRINT_BOTH("======================================\n");
    
    fprintf(stderr, "Counting k-mers ......\n");


    start = clock();
    if (count_bloom_shift) {
        h = count_files_filtered(argc - o.ind, argv + o.ind, k, p, block_size,
                                 n_thread, read_len, count_bloom_shift);
        if (!h) {
            fprintf(stderr, "ERROR: Bloom counting failed; check regular input files and memory budget\n");
            exit_status = 1;
            goto cleanup;
        }
    } else {
        fprintf(stderr, "Counting NGS file 1 ......\n");
        h = count_file(argv[o.ind], k, p, block_size, n_thread, read_len);
        if (!h) {
            fprintf(stderr, "ERROR: failed to load input file %s\n", argv[o.ind]);
            exit_status = 1;
            goto cleanup;
        }

        int f_num = argc - o.ind, c_f_n = 2;
        while(c_f_n <= f_num ) {
            fprintf(stderr, "Counting NGS file %d ......\n", c_f_n);
            if (!count_file2(argv[o.ind + c_f_n - 1], h, k, p, block_size, n_thread, read_len)) {
                fprintf(stderr, "ERROR: failed to load input file %s\n", argv[o.ind + c_f_n - 1]);
                exit_status = 1;
                goto cleanup;
            }
            c_f_n = c_f_n + 1;
        }

    }

    if (c4x_set_bloom(h, bloom_bits_per_key) < 0) {
        fprintf(stderr, "ERROR: failed to allocate Bloom query filter\n");
        exit_status = 1;
        goto cleanup;
    }
    PRINT_BOTH("Bloom query filter (-B): %d bits/key; %zu bitmap bytes\n",
               bloom_bits_per_key, h->bloom ? h->bloom->n_bytes : 0);

    end = clock();
    kc_dur = (double)(end - start);

    PRINT_BOTH("\nK-mer counting time: %f seconds\n",(kc_dur/CLOCKS_PER_SEC));
    if(mem_test){
            fprintf(stderr, "Memory test mode activated! Press ENTER to continue...\n");
            getchar();}
    fprintf(stderr, "Decoding strands from index %"PRIu64" to ",min_index);
    fprintf(stderr, "%"PRIu64" ......\n",max_index);

    uint64_t id, eid;
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
                if (codec) dbgps_codec_encode_index(codec, id, index_len, &eid);
                else eid = index_encode(id, index_len);

                c_index = find_paths_mapped(h, path, pr_uint64 | eid, mask, cov_cut, path_len, k, index_len, max_cov_ratio, skip_ratios, codec);
                if(c_index > -1) {
                    dec_strand_num++;
                    unsigned char DNA_seq[Max_Path_Len];

                    packed_dna_string(path[c_index], DNA_seq, path_len + index_len);
                    
                    if (output_prefix != NULL && strands_file != NULL) {
                        // Output to strands file if -o is specified
                        fprintf(strands_file, "%"PRIu64"\t%s\n", id, DNA_seq);
                    } else {
                        // Output to stdout if no output file is specified
                        printf("%"PRIu64"\t%s\n", id, DNA_seq);
                    }
                    free(path[c_index]);
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
                    if (codec) dbgps_codec_encode_index(codec, id, index_len, &eid);
                    else eid = index_encode(id, index_len);

                    c_index = find_paths_mapped(h, path, pr_uint64 | eid, mask, cov_cut, path_len, k, index_len, r, skip_ratios, codec);
                    if(c_index > -1) {
                        dec_strand_num++;
                        unsigned char DNA_seq[Max_Path_Len];

                        packed_dna_string(path[c_index], DNA_seq, path_len + index_len);
                        
                        if (output_prefix != NULL && strands_file != NULL) {
                            // Output to strands file if -o is specified
                            fprintf(strands_file, "%"PRIu64"\t%s\n", id, DNA_seq);
                        } else {
                            // Output to stdout if no output file is specified
                            printf("%"PRIu64"\t%s\n", id, DNA_seq);
                        }
                        free(path[c_index]);
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

cleanup:
    c4x_destroy(h);
    
    // Close output files if they were opened
    if (log_file != NULL) {
        fclose(log_file);
    }
    
    if (strands_file != NULL) {
        fclose(strands_file);
    }
    dbgps_codec_destroy(&mapping);
    
    #undef PRINT_BOTH
    return exit_status;
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
