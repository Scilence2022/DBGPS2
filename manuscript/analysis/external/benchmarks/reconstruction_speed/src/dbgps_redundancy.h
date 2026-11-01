#ifndef DBGPS_REDUNDANCY_H
#define DBGPS_REDUNDANCY_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Byte-level protection, independent of DNA mapping and graph traversal.
 * Configurations are immutable and can be shared by concurrent queries.
 * A zero-initialized configuration means none. */
typedef struct dbgps_redundancy_code dbgps_redundancy_code_t;
typedef struct {
    const dbgps_redundancy_code_t *code;
    unsigned parity; /* RS parity symbols (bytes); other codes ignore this field. */
} dbgps_redundancy_t;

#define DBGPS_REDUNDANCY_OPTIONS "none, crc16, crc32, crc32c, rs[:NB] (N bytes; default rs:16B; rs:N accepted), hamming"
#define DBGPS_REDUNDANCY_SPEC_SIZE 32

/* Parse strictly; rs:NB specifies bytes, rs:N remains a compatible shorthand,
 * and rs selects rs:16B. Lowercase b is rejected. Preserve config on error. */
int dbgps_redundancy_init(dbgps_redundancy_t *config, const char *spec,
                          char *error, size_t error_size);
const char *dbgps_redundancy_name(const dbgps_redundancy_t *config);
/* Canonical display specification, e.g. rs:8B. Uses caller-owned storage;
 * returns -1 if the buffer is missing/too small, otherwise 0. */
int dbgps_redundancy_format(const dbgps_redundancy_t *config, char *out, size_t out_size);
int dbgps_redundancy_enabled(const dbgps_redundancy_t *config);
int dbgps_redundancy_size(const dbgps_redundancy_t *config, size_t data_bytes,
                          size_t *parity_bytes);
/* Systematic wire format: data followed by parity. out may equal data. */
int dbgps_redundancy_encode(const dbgps_redundancy_t *config,
                            const unsigned char *data, size_t data_bytes,
                            unsigned char *out, size_t out_bytes);
/* Validate/correct a complete codeword in place. Return corrected symbol
 * count (bits for Hamming, bytes for RS), 0 if valid, -1 on failure.
 * The buffer is unchanged on failure. Beyond the correction radius, a code
 * can miscorrect to another codeword; no ECC guarantees detection there. */
int dbgps_redundancy_decode(const dbgps_redundancy_t *config,
                            unsigned char *word, size_t word_bytes, size_t data_bytes);

uint16_t dbgps_crc16_xmodem(const unsigned char *data, size_t bytes);

#ifdef __cplusplus
}
#endif
#endif
