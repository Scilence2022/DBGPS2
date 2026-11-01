#include "dbgps_redundancy.h"

#include <stdio.h>
#include <string.h>
#include <stdint.h>

/* To add a code, implement these three operations and add one registry row.
 * Layout, CLI, mapping, candidate validation and server queries use this same
 * dispatch, with no algorithm-specific branches in either graph engine. */
struct dbgps_redundancy_code {
    const char *name;
    int (*size)(const dbgps_redundancy_t *, size_t, size_t *);
    void (*encode)(const dbgps_redundancy_t *, unsigned char *, size_t, size_t);
    int (*decode)(const dbgps_redundancy_t *, unsigned char *, size_t, size_t);
    unsigned parameter; /* fixed checksum bytes; zero for parameterless ECC */
    uint32_t polynomial; /* reflected CRC32 polynomial */
};

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

uint16_t dbgps_crc16_xmodem(const unsigned char *data, size_t bytes)
{
    uint16_t crc = 0;
    for (size_t i = 0; i < bytes; ++i)
        crc = (uint16_t)((crc << 8) ^ crc16tab[(crc >> 8) ^ data[i]]);
    return crc;
}

static uint32_t crc32_value(const unsigned char *data, size_t bytes, uint32_t polynomial)
{
    uint32_t crc = UINT32_MAX;
    for (size_t i = 0; i < bytes; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit)
            crc = (crc >> 1) ^ ((crc & 1) ? polynomial : 0);
    }
    return crc ^ UINT32_MAX;
}

static int fixed_size(const dbgps_redundancy_t *c, size_t n, size_t *p)
{
    (void)n;
    *p = c->code->parameter;
    return 0;
}

static void checksum_encode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    uint32_t value = p == 2 ? dbgps_crc16_xmodem(w, n) :
        crc32_value(w, n, c->code->polynomial);
    for (size_t i = 0; i < p; ++i) w[n + i] = (unsigned char)(value >> (8 * (p - i - 1)));
}

static int checksum_decode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    unsigned char expected[4];
    uint32_t value = p == 2 ? dbgps_crc16_xmodem(w, n) :
        crc32_value(w, n, c->code->polynomial);
    for (size_t i = 0; i < p; ++i) expected[i] = (unsigned char)(value >> (8 * (p - i - 1)));
    return memcmp(w + n, expected, p) ? -1 : 0;
}

/* Shortened systematic RS over GF(256), primitive polynomial 0x11d,
 * alpha=2, generator roots alpha^0 .. alpha^(p-1), high coefficient first.
 * All tables and workspaces are local; no lazy global initialization. */
typedef struct { unsigned char exp[510], log[256]; } gf256_t;

static void gf_init(gf256_t *f)
{
    unsigned x = 1;
    memset(f, 0, sizeof(*f));
    for (unsigned i = 0; i < 255; ++i) {
        f->exp[i] = x;
        f->log[x] = i;
        x <<= 1;
        if (x & 256) x ^= 0x11d;
    }
    for (unsigned i = 255; i < 510; ++i) f->exp[i] = f->exp[i - 255];
}

static unsigned char gf_mul(const gf256_t *f, unsigned char a, unsigned char b)
{
    return a && b ? f->exp[f->log[a] + f->log[b]] : 0;
}

static unsigned char gf_div(const gf256_t *f, unsigned char a, unsigned char b)
{
    return a ? f->exp[f->log[a] + 255 - f->log[b]] : 0;
}

static int rs_size(const dbgps_redundancy_t *c, size_t n, size_t *p)
{
    if (c->parity < 1 || c->parity > 254 || !n || n > 255 - c->parity) return -1;
    *p = c->parity;
    return 0;
}

static void rs_encode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    (void)c;
    gf256_t f;
    gf_init(&f);
    unsigned char gen[256] = {1}, work[255];
    for (size_t i = 0; i < p; ++i) {
        gen[i + 1] = 0;
        for (size_t j = i + 1; j > 0; --j)
            gen[j] ^= gf_mul(&f, gen[j - 1], f.exp[i]);
    }
    memcpy(work, w, n);
    memset(work + n, 0, p);
    for (size_t i = 0; i < n; ++i) {
        unsigned char factor = work[i];
        for (size_t j = 1; j <= p; ++j) work[i + j] ^= gf_mul(&f, factor, gen[j]);
    }
    memcpy(w + n, work + n, p);
}

static int rs_syndromes(const gf256_t *f, const unsigned char *w, size_t n,
                         unsigned char *s, size_t p)
{
    unsigned nonzero = 0;
    for (size_t i = 0; i < p; ++i) {
        unsigned char value = 0;
        for (size_t j = 0; j < n; ++j) value = gf_mul(f, value, f->exp[i]) ^ w[j];
        nonzero |= s[i] = value;
    }
    return nonzero != 0;
}

static int rs_decode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    (void)c;
    gf256_t f;
    gf_init(&f);
    unsigned char s[254], locator[255] = {1}, previous[255] = {1};
    size_t length = n + p;
    if (!rs_syndromes(&f, w, length, s, p)) return 0;

    /* Berlekamp-Massey with ascending locator coefficients. */
    unsigned degree = 0, shift = 1;
    unsigned char last = 1;
    for (unsigned r = 0; r < p; ++r) {
        unsigned char delta = s[r];
        for (unsigned j = 1; j <= degree; ++j) delta ^= gf_mul(&f, locator[j], s[r - j]);
        if (!delta) { ++shift; continue; }
        unsigned char saved[255];
        memcpy(saved, locator, sizeof(saved));
        unsigned char scale = gf_div(&f, delta, last);
        for (unsigned j = 0; j + shift <= p; ++j)
            locator[j + shift] ^= gf_mul(&f, scale, previous[j]);
        if (2 * degree <= r) {
            degree = r + 1 - degree;
            memcpy(previous, saved, sizeof(previous));
            last = delta;
            shift = 1;
        } else ++shift;
    }
    if (!degree || 2 * degree > p) return -1;

    /* Chien search over positions actually present in the shortened word. */
    unsigned positions[127], found = 0;
    unsigned char locations[127];
    for (unsigned pos = 0; pos < length; ++pos) {
        unsigned power = (unsigned)length - 1 - pos;
        unsigned char inverse = f.exp[(255 - power) % 255];
        unsigned char value = locator[degree];
        for (unsigned j = degree; j > 0; --j) value = gf_mul(&f, value, inverse) ^ locator[j - 1];
        if (!value) {
            if (found >= degree) return -1;
            positions[found] = pos;
            locations[found++] = f.exp[power];
        }
    }
    if (found != degree) return -1;

    /* Error evaluator Omega=(S * locator) mod x^p; Forney for fcr=0.
     * Evaluate the formal derivative directly in characteristic two. */
    unsigned char omega[254] = {0}, repaired[255];
    for (unsigned i = 0; i < p; ++i)
        for (unsigned j = 0; j <= degree && i + j < p; ++j)
            omega[i + j] ^= gf_mul(&f, s[i], locator[j]);
    memcpy(repaired, w, length);
    for (unsigned i = 0; i < found; ++i) {
        unsigned char inverse = gf_div(&f, 1, locations[i]);
        unsigned char numerator = omega[p - 1], denominator = 0, power = 1;
        for (unsigned j = (unsigned)p - 1; j > 0; --j)
            numerator = gf_mul(&f, numerator, inverse) ^ omega[j - 1];
        for (unsigned j = 1; j <= degree; ++j) {
            if (j & 1) denominator ^= gf_mul(&f, locator[j], power);
            power = gf_mul(&f, power, inverse);
        }
        if (!denominator) return -1;
        repaired[positions[i]] ^= gf_mul(&f, locations[i], gf_div(&f, numerator, denominator));
    }
    if (rs_syndromes(&f, repaired, length, s, p)) return -1;
    memcpy(w, repaired, length);
    return (int)found;
}

/* Extended Hamming(13,8): data MSB first at positions 3,5,6,7,9,10,11,12;
 * even parity at 1,2,4,8 and overall even parity at 13. The five parity bits
 * per data byte are packed MSB first in an appended trailer. */
static const unsigned hamming_positions[8] = {3, 5, 6, 7, 9, 10, 11, 12};

static unsigned bit_parity(unsigned value)
{
    value ^= value >> 4;
    value ^= value >> 2;
    value ^= value >> 1;
    return value & 1;
}

static unsigned hamming_parity(unsigned char value)
{
    unsigned checks = 0;
    for (unsigned b = 0; b < 8; ++b)
        if (value & (1U << (7 - b))) checks ^= hamming_positions[b];
    return (checks << 1) | (bit_parity(value) ^ bit_parity(checks));
}

static int hamming_size(const dbgps_redundancy_t *c, size_t n, size_t *p)
{
    (void)c;
    if (n > (SIZE_MAX - 7) / 5) return -1;
    *p = (5 * n + 7) / 8;
    return 0;
}

static void hamming_encode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    (void)c;
    memset(w + n, 0, p);
    for (size_t i = 0; i < n; ++i) {
        unsigned checks = hamming_parity(w[i]);
        for (unsigned b = 0; b < 5; ++b) {
            size_t bit = 5 * i + b;
            w[n + bit / 8] |= ((checks >> (4 - b)) & 1) << (7 - bit % 8);
        }
    }
}

static int hamming_decode(const dbgps_redundancy_t *c, unsigned char *w, size_t n, size_t p)
{
    unsigned padding = (8 - (5 * n) % 8) % 8;
    if (p && (w[n + p - 1] & ((1U << padding) - 1))) return -1;
    int corrected = 0;
    /* Validate all codewords before changing any byte. */
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t i = 0; i < n; ++i) {
            unsigned checks = 0;
            for (unsigned b = 0; b < 5; ++b) {
                size_t bit = 5 * i + b;
                checks = (checks << 1) | ((w[n + bit / 8] >> (7 - bit % 8)) & 1);
            }
            unsigned syndrome = (checks ^ hamming_parity(w[i])) >> 1;
            unsigned odd = bit_parity(w[i]) ^ bit_parity(checks);
            if (syndrome && (!odd || syndrome > 12)) return -1;
            if (pass && odd) {
                ++corrected;
                for (unsigned b = 0; b < 8; ++b)
                    if (syndrome == hamming_positions[b]) w[i] ^= 1U << (7 - b);
            }
        }
    }
    hamming_encode(c, w, n, p);
    return corrected;
}

static const dbgps_redundancy_code_t registry[] = {
    {"crc16", fixed_size, checksum_encode, checksum_decode, 2, 0},
    {"crc32", fixed_size, checksum_encode, checksum_decode, 4, 0xedb88320U},
    {"crc32c", fixed_size, checksum_encode, checksum_decode, 4, 0x82f63b78U},
    {"rs", rs_size, rs_encode, rs_decode, 0, 0},
    {"hamming", hamming_size, hamming_encode, hamming_decode, 0, 0}
};

int dbgps_redundancy_init(dbgps_redundancy_t *config, const char *spec,
                          char *error, size_t error_size)
{
    dbgps_redundancy_t next = {0};
    if (!config || !spec) goto fail;
    if (!strcmp(spec, "none")) { *config = next; return 0; }
    for (size_t i = 0; i < sizeof(registry) / sizeof(*registry); ++i) {
        size_t len = strlen(registry[i].name);
        if (strncmp(spec, registry[i].name, len)) continue;
        next.code = &registry[i];
        if (!strcmp(registry[i].name, "rs")) {
            next.parity = 16;
            if (spec[len]) {
                if (spec[len++] != ':' || !spec[len]) goto fail;
                unsigned value = 0;
                for (; spec[len] >= '0' && spec[len] <= '9'; ++len) {
                    if (value > 254) goto fail;
                    value = value * 10 + (spec[len] - '0');
                }
                if (spec[len] == 'B') ++len;
                if (spec[len] || value < 1 || value > 254) goto fail;
                next.parity = value;
            }
        } else if (spec[len]) continue;
        *config = next;
        return 0;
    }
fail:
    if (error && error_size) snprintf(error, error_size, "invalid redundancy code; use " DBGPS_REDUNDANCY_OPTIONS "; RS N must be 1..254 bytes (uppercase B)");
    return -1;
}

const char *dbgps_redundancy_name(const dbgps_redundancy_t *c)
{
    return c && c->code ? c->code->name : "none";
}

int dbgps_redundancy_format(const dbgps_redundancy_t *c, char *out, size_t out_size)
{
    if (!out || !out_size) return -1;
    const char *name = dbgps_redundancy_name(c);
    int length = !strcmp(name, "rs") ? snprintf(out, out_size, "rs:%uB", c->parity) :
        snprintf(out, out_size, "%s", name);
    return length < 0 || (size_t)length >= out_size ? -1 : 0;
}

int dbgps_redundancy_enabled(const dbgps_redundancy_t *c)
{
    return c && c->code;
}

int dbgps_redundancy_size(const dbgps_redundancy_t *c, size_t n, size_t *p)
{
    if (!p) return -1;
    if (!dbgps_redundancy_enabled(c)) *p = 0;
    else if (c->code->size(c, n, p) < 0) return -1;
    return n <= SIZE_MAX - *p ? 0 : -1;
}

int dbgps_redundancy_encode(const dbgps_redundancy_t *c, const unsigned char *data,
                            size_t n, unsigned char *out, size_t out_bytes)
{
    size_t p;
    if (!out || (!data && n) || dbgps_redundancy_size(c, n, &p) < 0 || out_bytes != n + p) return -1;
    if (n) memmove(out, data, n);
    if (dbgps_redundancy_enabled(c)) c->code->encode(c, out, n, p);
    return 0;
}

int dbgps_redundancy_decode(const dbgps_redundancy_t *c, unsigned char *w,
                            size_t word_bytes, size_t n)
{
    size_t p;
    if (!w || dbgps_redundancy_size(c, n, &p) < 0 || word_bytes != n + p) return -1;
    return dbgps_redundancy_enabled(c) ? c->code->decode(c, w, n, p) : 0;
}
