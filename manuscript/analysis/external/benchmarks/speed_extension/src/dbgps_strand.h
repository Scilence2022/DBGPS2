#ifndef DBGPS_STRAND_H
#define DBGPS_STRAND_H

#include "dbgps_codec.h"
#include "dbgps_redundancy.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    unsigned index_bases, index_bytes, data_bytes, code_bytes;
    unsigned parity_bytes, parity_bases, padding_bits, path_bases, total_bases;
} dbgps_strand_layout_t;

/* Compose byte-level redundancy with fixed-block DNA mapping. Lengths are in
 * bases and describe data only; parity and mapping padding are derived. */
int dbgps_strand_layout(const dbgps_codec_t *codec, const dbgps_redundancy_t *redundancy,
                        int index_bases, int payload_bases, dbgps_strand_layout_t *out);
/* Decode, validate/correct, and re-encode the candidate, preserving its seeded
 * index. Distinct graph paths can yield the same repaired DNA. Buffers must be
 * separate: input/output need ceil(total_bases/4) bytes; scratch code_bytes.
 * A forbidden mapping word or nonzero padding is rejected before ECC. */
int dbgps_strand_check(const dbgps_codec_t *codec, const dbgps_redundancy_t *redundancy,
                       const dbgps_strand_layout_t *layout, const unsigned char *candidate,
                       unsigned char *corrected, unsigned char *scratch);

#ifdef __cplusplus
}
#endif
#endif
