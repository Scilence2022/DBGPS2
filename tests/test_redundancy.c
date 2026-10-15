#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "dbgps_redundancy.h"
#include "dbgps_strand.h"

static dbgps_redundancy_t config(const char *spec)
{
    dbgps_redundancy_t c = {0};
    char error[256];
    assert(dbgps_redundancy_init(&c, spec, error, sizeof(error)) == 0);
    return c;
}

static void test_checksums(void)
{
    const char *names[] = {"crc16", "crc32", "crc32c"};
    const unsigned char expected[][4] = {{0x31, 0xc3}, {0xcb, 0xf4, 0x39, 0x26}, {0xe3, 0x06, 0x92, 0x83}};
    for (unsigned i = 0; i < 3; ++i) {
        dbgps_redundancy_t c = config(names[i]);
        size_t p;
        assert(dbgps_redundancy_size(&c, 9, &p) == 0);
        unsigned char word[13], saved[13];
        assert(dbgps_redundancy_encode(&c, (const unsigned char *)"123456789", 9, word, 9 + p) == 0);
        assert(!memcmp(word + 9, expected[i], p));
        assert(dbgps_redundancy_decode(&c, word, 9 + p, 9) == 0);
        for (size_t j = 0; j < 9 + p; ++j) {
            word[j] ^= 1;
            memcpy(saved, word, 9 + p);
            assert(dbgps_redundancy_decode(&c, word, 9 + p, 9) == -1);
            assert(!memcmp(saved, word, 9 + p));
            word[j] ^= 1;
        }
        assert(dbgps_redundancy_size(&c, 0, &p) == 0);
        assert(dbgps_redundancy_encode(&c, NULL, 0, word, p) == 0);
        assert(dbgps_redundancy_decode(&c, word, p, 0) == 0);
    }
}

static unsigned random_state = 192193;
static unsigned next_random(void)
{
    random_state = random_state * 1664525U + 1013904223U;
    return random_state;
}

static void test_rs(void)
{
    /* Independent reedsolo RSCodec(10, prim=0x11d, generator=2, fcr=0). */
    const unsigned char vector[] = {0x68,0x65,0x6c,0x6c,0x6f,0x20,0x77,0x6f,0x72,0x6c,0x64,
                                    0xed,0x25,0x54,0xc4,0xfd,0xfd,0x89,0xf3,0xa8,0xaa};
    unsigned char word[255], damaged[255];
    dbgps_redundancy_t c = config("rs:10B");
    assert(dbgps_redundancy_encode(&c, vector, 11, word, sizeof(vector)) == 0);
    assert(!memcmp(word, vector, sizeof(vector)));
    for (unsigned parity = 1; parity <= 254; ++parity) {
        char spec[32];
        snprintf(spec, sizeof(spec), "rs:%u", parity);
        dbgps_redundancy_t legacy = config(spec);
        snprintf(spec, sizeof(spec), "rs:%uB", parity);
        c = config(spec);
        assert(c.code == legacy.code && c.parity == legacy.parity);
        char label[DBGPS_REDUNDANCY_SPEC_SIZE];
        assert(dbgps_redundancy_format(&c, label, sizeof(label)) == 0);
        assert(!strcmp(label, spec));
        for (unsigned trial = 0; trial < 3; ++trial) {
            size_t n = trial == 0 ? 1 : trial == 1 ? 255 - parity : 1 + next_random() % (255 - parity);
            size_t length = n + parity;
            for (size_t i = 0; i < n; ++i) word[i] = next_random() >> 24;
            assert(dbgps_redundancy_encode(&c, word, n, word, length) == 0);
            assert(dbgps_redundancy_decode(&c, word, length, n) == 0);
            memcpy(damaged, word, length);
            unsigned positions[255];
            for (unsigned i = 0; i < length; ++i) positions[i] = i;
            unsigned errors = parity / 2;
            for (unsigned i = 0; i < errors; ++i) {
                unsigned pick = i + next_random() % (length - i);
                unsigned at = positions[pick];
                positions[pick] = positions[i];
                damaged[at] ^= 1 + next_random() % 255;
            }
            assert(dbgps_redundancy_decode(&c, damaged, length, n) == (int)errors);
            assert(!memcmp(damaged, word, length));
        }
    }
    c = config("rs:4");
    memset(word, 0x53, 20);
    assert(dbgps_redundancy_encode(&c, word, 20, word, 24) == 0);
    word[0] ^= 1; word[1] ^= 2; word[2] ^= 4;
    memcpy(damaged, word, 24);
    assert(dbgps_redundancy_decode(&c, word, 24, 20) == -1);
    assert(!memcmp(damaged, word, 24));
    size_t p;
    assert(dbgps_redundancy_size(&c, 252, &p) == -1);
    assert(dbgps_redundancy_size(&c, 0, &p) == -1);
}

static void flip(unsigned char *w, unsigned bit)
{
    w[bit / 8] ^= 1U << (7 - bit % 8);
}

static void test_hamming(void)
{
    dbgps_redundancy_t c = config("hamming");
    /* All values, all thirteen single-bit errors, all double-bit errors. */
    for (unsigned value = 0; value < 256; ++value) {
        unsigned char original[2] = {value, 0}, damaged[2], saved[2];
        assert(dbgps_redundancy_encode(&c, original, 1, original, 2) == 0);
        assert(dbgps_redundancy_decode(&c, original, 2, 1) == 0);
        for (unsigned a = 0; a < 13; ++a) {
            memcpy(damaged, original, 2);
            flip(damaged, a);
            assert(dbgps_redundancy_decode(&c, damaged, 2, 1) == 1);
            assert(!memcmp(damaged, original, 2));
            for (unsigned b = a + 1; b < 13; ++b) {
                memcpy(damaged, original, 2);
                flip(damaged, a); flip(damaged, b);
                memcpy(saved, damaged, 2);
                assert(dbgps_redundancy_decode(&c, damaged, 2, 1) == -1);
                assert(!memcmp(damaged, saved, 2));
            }
        }
        original[1] |= 1; /* unused trailer bits cannot carry data */
        assert(dbgps_redundancy_decode(&c, original, 2, 1) == -1);
    }
    unsigned char w[13], saved[13];
    for (unsigned i = 0; i < 8; ++i) w[i] = i * 31;
    assert(dbgps_redundancy_encode(&c, w, 8, w, 13) == 0);
    memcpy(saved, w, 13);
    for (unsigned i = 0; i < 8; ++i) w[i] ^= 0x80;
    assert(dbgps_redundancy_decode(&c, w, 13, 8) == 8);
    assert(!memcmp(w, saved, 13));
    w[0] ^= 0x80; w[7] ^= 3;
    memcpy(saved, w, 13);
    assert(dbgps_redundancy_decode(&c, w, 13, 8) == -1);
    assert(!memcmp(w, saved, 13));
}

static void test_configuration(void)
{
    const char *invalid[] = {"", "bch", "rs:", "rs:0", "rs:255", "rs:-1", "rs:+4",
        "rs:4x", "rs:99999999999999999999999", "crc16:2", "none:0", "hamming:4",
        "rs:8b", "rs:8bits", "rs:B", "rs:0B", "rs:255B", "rs:-8B", "rs:+8B",
        "rs:8BB", "rs:8Bjunk", "rs:8 B", "rs:8B ", "rs:99999999999999999999999B"};
    dbgps_redundancy_t c = config("rs");
    assert(c.parity == 16);
    char label[DBGPS_REDUNDANCY_SPEC_SIZE];
    assert(dbgps_redundancy_format(&c, label, sizeof(label)) == 0);
    assert(!strcmp(label, "rs:16B"));
    char tiny[3];
    assert(dbgps_redundancy_format(&c, tiny, sizeof(tiny)) == -1);
    assert(dbgps_redundancy_format(&c, NULL, 0) == -1);
    for (unsigned i = 0; i < sizeof(invalid) / sizeof(*invalid); ++i) {
        char error[256];
        assert(dbgps_redundancy_init(&c, invalid[i], error, sizeof(error)) == -1);
        assert(strstr(error, "invalid"));
        assert(c.parity == 16 && !strcmp(dbgps_redundancy_name(&c), "rs"));
    }
    c = config("none");
    assert(dbgps_redundancy_format(&c, label, sizeof(label)) == 0 && !strcmp(label, "none"));
    c = config("crc16");
    assert(dbgps_redundancy_format(&c, label, sizeof(label)) == 0 && !strcmp(label, "crc16"));
    c = config("none");
    size_t p;
    assert(!dbgps_redundancy_enabled(&c));
    assert(dbgps_redundancy_size(&c, SIZE_MAX, &p) == 0 && p == 0);
    c = config("hamming");
    assert(dbgps_redundancy_size(&c, SIZE_MAX, &p) == -1);
    c = config("crc16");
    assert(dbgps_redundancy_size(&c, SIZE_MAX, &p) == -1);
    unsigned char w[8] = {0};
    assert(dbgps_redundancy_encode(&c, w, 4, w, 5) == -1);
    assert(dbgps_redundancy_decode(&c, w, 5, 4) == -1);
}

static void test_index_guard(void)
{
    dbgps_redundancy_t c = config("rs:4");
    dbgps_strand_layout_t layout;
    assert(dbgps_strand_layout(NULL, &c, 4, 12, &layout) == 0);
    unsigned char raw[8] = {7, 25, 81, 189}, dna[8], corrected[8], scratch[8];
    assert(dbgps_redundancy_encode(&c, raw, 4, raw, 8) == 0);
    assert(dbgps_codec_encode(NULL, raw, 8, dna, 32) == 0);
    dna[0] ^= 1; /* correction must not silently change the queried index */
    assert(dbgps_strand_check(NULL, &c, &layout, dna, corrected, scratch) == -1);
    dna[0] ^= 1; dna[2] ^= 1;
    assert(dbgps_strand_check(NULL, &c, &layout, dna, corrected, scratch) == 1);
    assert(!memcmp(corrected, raw, 8));
}

int main(void)
{
    test_configuration();
    test_checksums();
    test_rs();
    test_hamming();
    test_index_guard();
    puts("redundancy: checksum vectors, RS1..254, correction limits, exhaustive SECDED, atomic failures and index guard passed");
    return 0;
}
