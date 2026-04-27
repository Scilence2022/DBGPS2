#ifndef DBGPS3_CORE_H
#define DBGPS3_CORE_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Constants */
#define DBGPS3_KC_BITS 14
#define DBGPS3_KC_MAX ((1<<DBGPS3_KC_BITS) - 1)
#define DBGPS3_MAX_PATH_NUM 1050
#define DBGPS3_MAX_PATH_LEN 300

/* Translation tables */
extern const unsigned char dbgps3_seq_nt4_table[256];
extern const unsigned char dbgps3_nt4_seq_table[4];

/* Core graph object */
typedef struct dbgps3_graph_s {
    void *_h;          /* internal hash table (kc_c4x_t*), opaque */
    int k;             /* k-mer size */
    int p;             /* prefix bits */
    uint64_t mask;     /* k-mer bit mask */
    int index_len;     /* index length in bases */
    int index_byte_num;/* index length in bytes */
    int path_len;      /* payload path length */
    int n_thread;      /* number of threads */
    int block_size;    /* pipeline block size */
    int read_len;      /* read length for counting */
} dbgps3_graph_t;

/* ---- Graph lifecycle ---- */

/**
 * Load a de Bruijn graph from one or more FASTQ/FASTA files.
 * Returns a new graph object, or NULL on failure.
 */
dbgps3_graph_t *dbgps3_load(int k, int p, int n_thread, int block_size,
                             int read_len, int n_files, char **files);

/**
 * Destroy a graph object and free all associated memory.
 */
void dbgps3_destroy(dbgps3_graph_t *g);

/**
 * Load additional FASTQ/FASTA files into an existing graph.
 * Returns 0 on success, -1 on failure.
 */
int dbgps3_load_more(dbgps3_graph_t *g, int n_files, char **files);

/* ---- Query API ---- */

/**
 * Query the coverage of a k-mer (given as a min-hash key).
 * Returns 0 if the k-mer is not in the graph.
 */
unsigned short dbgps3_kmer_cov(const dbgps3_graph_t *g, uint64_t kmer);

/**
 * Find a path in the de Bruijn graph starting from start_kmer.
 * On success, writes the decoded DNA sequence into result_buf and returns 0.
 * On failure (no unique path found), returns -1 and result_buf is unchanged.
 *
 * @param g            Graph object
 * @param start_kmer   Starting k-mer (uint64 encoded)
 * @param cut_cov      Minimum coverage cutoff
 * @param max_cov_ratio Maximum coverage ratio for branching
 * @param skip_ratios  Number of initial positions to skip ratio checks
 * @param path_len     Payload path length (in bases)
 * @param index_byte_num Number of bytes in the index encoding
 * @param result_buf   Output buffer for the decoded DNA sequence
 * @param result_buf_len Size of result_buf (must be >= path_len + index_byte_num*4 + 1)
 */
int dbgps3_find_path(const dbgps3_graph_t *g, uint64_t start_kmer,
                     int cut_cov, double max_cov_ratio, int skip_ratios,
                     int path_len, int index_byte_num,
                     char *result_buf, int result_buf_len);

/* ---- Utility functions ---- */

/** Encode an ACGT sequence to a uint64 k-mer (forward strand only). */
uint64_t dbgps3_actgkmer_uint64(const char *seq, int km_len);

/** Encode an ACGT sequence to a uint64 k-mer (min-hash of forward/reverse). */
uint64_t dbgps3_actgkmer_hashkey(const char *seq, int km_len);

/** Compute the primer k-mer as uint64 from a primer sequence. */
uint64_t dbgps3_primer_uint64(const char *pr, int pr_len, int index_len, int km_len);

/** Encode a decimal index into DNA encoding. */
uint64_t dbgps3_index_encode(uint64_t index, unsigned int index_len);

/** Compute reverse complement of a k-mer. */
uint64_t dbgps3_comp_rev(uint64_t x, unsigned char km_len);

/** Compute min-hash key (min of forward and reverse complement). */
uint64_t dbgps3_min_hash_key(uint64_t x, unsigned char km_len);

/** Decode a uint64 k-mer to an ACGT sequence string. */
void dbgps3_uint64_acgt(uint64_t key, char *seq, int km_len);

/* ---- Graph accessors ---- */

int dbgps3_get_k(const dbgps3_graph_t *g);
uint64_t dbgps3_get_mask(const dbgps3_graph_t *g);
int dbgps3_get_index_len(const dbgps3_graph_t *g);
int dbgps3_get_index_byte_num(const dbgps3_graph_t *g);
int dbgps3_get_path_len(const dbgps3_graph_t *g);
int dbgps3_get_p(const dbgps3_graph_t *g);
long dbgps3_get_kmer_count(const dbgps3_graph_t *g);

#ifdef __cplusplus
}
#endif

#endif /* DBGPS3_CORE_H */
