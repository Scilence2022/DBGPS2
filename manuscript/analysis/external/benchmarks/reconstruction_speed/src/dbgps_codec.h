#ifndef DBGPS_CODEC_H
#define DBGPS_CODEC_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Immutable after loading; graph keys always retain A=0,C=1,G=2,T=3. */
typedef struct {
    unsigned bases, bits;
    uint16_t *encode;
    uint32_t *decode; /* UINT32_MAX denotes a forbidden DNA word. */
    int byte_fast;
    unsigned char encode_byte[256], decode_byte[256];
} dbgps_codec_t;

/* Load into an empty codec. On error, leave it unchanged and describe the
   failure in error. Destroy a loaded codec before loading another codebook. */
int dbgps_codec_load(dbgps_codec_t *codec, const char *filename, char *error, size_t error_size);
void dbgps_codec_destroy(dbgps_codec_t *codec);
int dbgps_codec_encode_index(const dbgps_codec_t *codec, uint64_t index,
                             unsigned index_bases, uint64_t *encoded);
/* DNA is packed MSB first, four canonical bases per byte. Final unused slots
   are zero. Decode exactly out_bytes and require zero trailing padding bits.
   A NULL codec selects canonical A=00,C=01,G=10,T=11 packing. */
int dbgps_codec_decode(const dbgps_codec_t *codec, const unsigned char *dna,
                       unsigned dna_bases, unsigned char *out, unsigned out_bytes);
/* Encode complete bytes, appending zero bits to fill the final mapping block.
   dna_bases must be the exact resulting length; NULL selects canonical packing. */
int dbgps_codec_encode(const dbgps_codec_t *codec, const unsigned char *raw,
                       unsigned raw_bytes, unsigned char *dna, unsigned dna_bases);

#ifdef __cplusplus
}
#endif

#endif
