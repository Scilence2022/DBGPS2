#ifndef DBGPS_CORE_H
#define DBGPS_CORE_H

#include <stdint.h>
#include <stddef.h>
#include "dbgps_codec.h"
#include "dbgps_strand.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Constants */
#define DBGPS_KC_BITS 14
#define DBGPS_KC_MAX ((1<<DBGPS_KC_BITS) - 1)
#define DBGPS_MAX_PATH_NUM 1050
#define DBGPS_MAX_PATH_LEN 300
#define DBGPS_DEFAULT_REDUNDANCY "crc16"

/* Translation tables */
extern const unsigned char dbgps_seq_nt4_table[256];
extern const unsigned char dbgps_nt4_seq_table[4];

/* Core graph object */
typedef struct dbgps_graph_s {
    void *_h;          /* internal hash table (kc_c4x_t*), opaque */
    int k;             /* k-mer size */
    int p;             /* prefix bits */
    uint64_t mask;     /* k-mer bit mask */
    int index_len;     /* index length in bases */
    int index_byte_num;/* index length in decoded data bytes */
    int path_len;      /* payload path length */
    int n_thread;      /* number of threads */
    int block_size;    /* pipeline block size */
    int read_len;      /* read length for counting */
    int bloom_bits_per_key; /* optional query filter budget; 0 disables */
    int count_bloom_shift; /* two-pass singleton filtering; 0 disables */
    char **_source_files;  /* owned replay list for filtered graph/load */
    int _n_source_files;
    dbgps_codec_t codec; /* owned; bases=0 means direct canonical packing */
    dbgps_redundancy_t redundancy; /* immutable query configuration; zero=none */
} dbgps_graph_t;

/* ---- Graph lifecycle ---- */

/**
 * Load a de Bruijn graph from one or more FASTQ/FASTA files.
 * Returns a new graph object, or NULL on failure. Loaded graphs default to
 * CRC16 protection; call dbgps_set_redundancy(g, "none", ...) to disable it.
 */
dbgps_graph_t *dbgps_load(int k, int p, int n_thread, int block_size,
                             int read_len, int n_files, char **files);

/** Optional two-pass counting: count_bloom_shift=0 preserves all k-mers;
 * 23..40 allocates 2^shift bitmap bits for candidate selection, then counts
 * candidates exactly and removes singletons. Inputs must be regular files
 * available for replay, including subsequent load_more calls. */
dbgps_graph_t *dbgps_load_filtered(int k, int p, int n_thread, int block_size,
                                     int read_len, int n_files, char **files,
                                     int count_bloom_shift);

/**
 * Destroy a graph object and free all associated memory.
 */
void dbgps_destroy(dbgps_graph_t *g);

/** Configure an optional exact-query prefilter (0 disables; 1..64 bits/key).
 * Build from all existing k-mers; preserve the previous filter on failure.
 * Call only while queries are stopped. No coverage counts are changed. */
int dbgps_set_bloom(dbgps_graph_t *g, int bits_per_key);

/** Bitmap bytes currently allocated (0 when disabled). */
size_t dbgps_bloom_bytes(const dbgps_graph_t *g);

/** Load an optional codebook (NULL resets canonical encoding). Set compatible
 * index_len/path_len defaults before calling; index_byte_num is recalculated.
 * On error, preserve the previous codec. Configure only before starting queries;
 * the graph owns the codebook until destruction. */
int dbgps_set_mapping(dbgps_graph_t *g, const char *filename,
                       char *error, size_t error_size);

/** Select byte-level redundancy before starting queries. Validate the default
 * layout and preserve the old configuration on failure. NULL spec is invalid. */
int dbgps_set_redundancy(dbgps_graph_t *g, const char *spec,
                          char *error, size_t error_size);

/** Validate index/payload lengths in bases and return decoded index bytes. */
int dbgps_validate_layout(const dbgps_graph_t *g, int index_len, int path_len,
                           int *index_bytes);

/** Encode an index with the graph's selected codec; reject capacity overflow. */
int dbgps_encode_index(const dbgps_graph_t *g, uint64_t index, int index_len,
                        uint64_t *encoded);

/** Parse an exact unsigned decimal index (also used by CLI and JSON strings). */
int dbgps_parse_index(const char *text, uint64_t *index);

/**
 * Load additional FASTQ/FASTA files into an existing graph.
 * Returns 0 on success, -1 on failure.
 * Queries must be stopped. Rebuild an enabled Bloom filter after loading;
 * on failure, partial counts may remain and filtering is disabled so queries
 * continue to use exact counts without a stale filter.
 * In two-pass counting mode, rebuild from all recorded inputs plus new files
 * and replace the graph only on success; preserve the old graph on failure.
 */
int dbgps_load_more(dbgps_graph_t *g, int n_files, char **files);

/* ---- Query API ---- */

/**
 * Query the coverage of a k-mer (given as a min-hash key).
 * Returns 0 if the k-mer is not in the graph.
 */
unsigned short dbgps_kmer_cov(const dbgps_graph_t *g, uint64_t kmer);

/**
 * Find a path in the de Bruijn graph starting from start_kmer.
 * On success, writes the decoded DNA sequence into result_buf and returns 0.
 * On failure (no unique validated/corrected result, workspace limit, invalid parameters or allocation
 * failure), returns -1 and result_buf is unchanged.
 *
 * @param g            Graph object
 * @param start_kmer   Starting k-mer (uint64 encoded)
 * @param cut_cov      Minimum coverage cutoff
 * @param max_cov_ratio Maximum coverage ratio for branching (finite and >= 1)
 * @param skip_ratios  Number of initial positions to skip ratio checks
 * @param path_len     Payload length in bases, from 0 to DBGPS_MAX_PATH_LEN;
 *                     must contain complete codewords and decode to whole bytes.
 *                     Selected redundancy is appended automatically.
 * @param index_byte_num Number of decoded data bytes in the index
 * @param result_buf   Output buffer for the decoded DNA sequence
 * @param result_buf_len Size of result_buf (>= index + payload + redundancy bases + 1)
 */
int dbgps_find_path(const dbgps_graph_t *g, uint64_t start_kmer,
                     int cut_cov, double max_cov_ratio, int skip_ratios,
                     int path_len, int index_byte_num,
                     char *result_buf, int result_buf_len);

/* ---- Utility functions ---- */

/** Encode an ACGT sequence to a uint64 k-mer (forward strand only). */
uint64_t dbgps_actgkmer_uint64(const char *seq, int km_len);

/** Encode an ACGT sequence to a uint64 k-mer (min-hash of forward/reverse). */
uint64_t dbgps_actgkmer_hashkey(const char *seq, int km_len);

/** Compute the primer k-mer as uint64 from a primer sequence. */
uint64_t dbgps_primer_uint64(const char *pr, int pr_len, int index_len, int km_len);

/** Use index bits directly as canonical DNA (A=00,C=01,G=10,T=11).
 * Use dbgps_encode_index() when a graph-specific codebook is needed. */
uint64_t dbgps_index_encode(uint64_t index, unsigned int index_len);

/** Compute reverse complement of a k-mer. */
uint64_t dbgps_comp_rev(uint64_t x, unsigned char km_len);

/** Compute min-hash key (min of forward and reverse complement). */
uint64_t dbgps_min_hash_key(uint64_t x, unsigned char km_len);

/** Decode a uint64 k-mer to an ACGT sequence string. */
void dbgps_uint64_acgt(uint64_t key, char *seq, int km_len);

/* ---- Graph accessors ---- */

int dbgps_get_k(const dbgps_graph_t *g);
uint64_t dbgps_get_mask(const dbgps_graph_t *g);
int dbgps_get_index_len(const dbgps_graph_t *g);
int dbgps_get_index_byte_num(const dbgps_graph_t *g);
int dbgps_get_path_len(const dbgps_graph_t *g);
int dbgps_get_p(const dbgps_graph_t *g);
long dbgps_get_kmer_count(const dbgps_graph_t *g);

#ifdef __cplusplus
}
#endif

#endif /* DBGPS_CORE_H */
