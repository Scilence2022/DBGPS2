#include "dbgps_strand.h"

#include <limits.h>
#include <string.h>

int dbgps_strand_layout(const dbgps_codec_t *codec, const dbgps_redundancy_t *redundancy,
                        int index_bases, int payload_bases, dbgps_strand_layout_t *out)
{
    unsigned bases = codec ? codec->bases : 1, bits = codec ? codec->bits : 2;
    if (!out || !bases || !bits || index_bases < 0 || payload_bases < 0 ||
        index_bases % bases || payload_bases % bases) return -1;
    uint64_t index_bits = (uint64_t)index_bases / bases * bits;
    uint64_t payload_bits = (uint64_t)payload_bases / bases * bits;
    if (index_bits % 8 || payload_bits % 8 || (index_bits + payload_bits) / 8 > UINT_MAX) return -1;
    size_t data_bytes = (index_bits + payload_bits) / 8, parity;
    if (dbgps_redundancy_size(redundancy, data_bytes, &parity) < 0 ||
        parity > UINT_MAX / 8 || data_bytes > UINT_MAX - parity) return -1;
    uint64_t parity_bases = ((uint64_t)parity * 8 + bits - 1) / bits * bases;
    uint64_t total = (uint64_t)index_bases + payload_bases + parity_bases;
    if (total > INT_MAX) return -1;
    dbgps_strand_layout_t next = {
        (unsigned)index_bases, (unsigned)(index_bits / 8), (unsigned)data_bytes,
        (unsigned)(data_bytes + parity), (unsigned)parity, (unsigned)parity_bases,
        (unsigned)((bits - parity * 8 % bits) % bits),
        (unsigned)(payload_bases + parity_bases), (unsigned)total
    };
    *out = next;
    return 0;
}

int dbgps_strand_check(const dbgps_codec_t *codec, const dbgps_redundancy_t *redundancy,
                       const dbgps_strand_layout_t *layout, const unsigned char *candidate,
                       unsigned char *corrected, unsigned char *scratch)
{
    if (dbgps_codec_decode(codec, candidate, layout->total_bases, scratch, layout->code_bytes) < 0)
        return -1;
    int repaired = dbgps_redundancy_decode(redundancy, scratch, layout->code_bytes, layout->data_bytes);
    if (repaired < 0) return -1;
    if (!repaired) {
        memcpy(corrected, candidate, (layout->total_bases + 3) / 4);
        return 0;
    }
    if (dbgps_codec_encode(codec, scratch, layout->code_bytes,
                            corrected, layout->total_bases) < 0) return -1;
    unsigned bytes = layout->index_bases / 4, remainder = layout->index_bases % 4;
    if (memcmp(candidate, corrected, bytes) || (remainder &&
        ((candidate[bytes] ^ corrected[bytes]) & (0xffU << (8 - 2 * remainder))))) return -1;
    return repaired;
}
