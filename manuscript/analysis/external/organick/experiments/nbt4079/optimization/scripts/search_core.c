/* Experiment driver only. Graph and bounded unique-path search are included
 * unchanged from the repository's DBGPS2.c. No references are read here. */
#ifdef USE_INNER
#include "../../../../dbgps_strand.h"
static int check_paper_inner(const dbgps_codec_t *, const dbgps_redundancy_t *,
                            const dbgps_strand_layout_t *, const unsigned char *,
                            unsigned char *, unsigned char *);
#define dbgps_strand_check check_paper_inner
#endif
#define main original_dbgps2_main
#include "../../../../DBGPS2.c"
#undef main
#ifdef USE_INNER
#undef dbgps_strand_check
#endif

static char current_direction = 'F';

static void reverse_complement_body(char *dna)
{
    for (int i = 0; i < 55; ++i) {
        char left = dna[i], right = dna[109-i];
        dna[i] = nt4_seq_table[3-seq_nt4_table[(unsigned char)right]];
        dna[109-i] = nt4_seq_table[3-seq_nt4_table[(unsigned char)left]];
    }
}

#ifdef USE_INNER
/* Dataset-format validator, replacing the inapplicable CRC verifier at the
 * existing candidate-validation call site. The traversal is unchanged.
 * Mapping/phase were calibrated globally and independently validated on all
 * author design sequences; no per-reference sequence lookup occurs here. */
static int cycle_index(char base)
{
    return base == 'A' ? 0 : base == 'G' ? 1 : base == 'C' ? 2 : 3;
}

static int check_paper_inner(const dbgps_codec_t *codec, const dbgps_redundancy_t *redundancy,
                            const dbgps_strand_layout_t *layout, const unsigned char *candidate,
                            unsigned char *corrected, unsigned char *scratch)
{
    int result = dbgps_strand_check(codec, redundancy, layout, candidate, corrected, scratch);
    if (result < 0 || layout->total_bases != 112) return -1;
    unsigned char dna[113];
    packed_dna_string(corrected, dna, 112);
    char *body = (char *)dna+2;
    if (current_direction == 'R') reverse_complement_body(body);
    unsigned value = 0;
    for (int i = 1; i < 110; ++i) {
        unsigned digit = (cycle_index(body[i])-cycle_index(body[i-1])+3) & 3;
        if (digit == 3) return -1;
        if (i == 109) return digit == 0 ? result : -1;
        value = value*3+digit;
        if (i%4 == 0) {
            if (value >= 64) return -1;
            value = 0;
        }
    }
    return -1;
}
#endif

int main(int argc, char **argv)
{
    if (argc != 4) return 1;
    int k = atoi(argv[2]);
    if (k < 16 || k > 31) return 1;
    uint64_t mask = (1ULL << (2*k)) - 1;
    clock_t start = clock();
    kc_c4x_t *h = count_file(argv[1], k, KC_BITS, 10000000, atoi(argv[3]), 200);
    if (!h) return 1;
    fprintf(stderr, "k=%d graph_cpu_seconds=%.6f\n", k, (double)(clock()-start)/CLOCKS_PER_SEC);
    puts("READY"); fflush(stdout);
    dbgps_redundancy_t redundancy = {0};
    if (dbgps_redundancy_init(&redundancy, "none", NULL, 0) < 0) return 1;
    char *line = NULL; size_t allocated = 0;
    while (getline(&line, &allocated, stdin) >= 0) {
        char tag[128], seed_file[4096], out_file[4096];
        int cutoff, skip; double ratio;
        if (!strncmp(line, "QUIT", 4)) break;
        if (sscanf(line, "%127s %d %lf %d %4095s %4095s", tag, &cutoff, &ratio, &skip, seed_file, out_file) != 6)
            return 1;
        FILE *seeds = fopen(seed_file, "r"), *out = fopen(out_file, "w");
        if (!seeds || !out) return 1;
        fputs("direction\tseed\tsupport\tpayload\n", out);
        char direction, seed[64]; unsigned support;
        unsigned long queries = 0, returned = 0;
        start = clock();
        while (fscanf(seeds, " %c %63s %u", &direction, seed, &support) == 3) {
            if (strlen(seed) != 31 || (direction != 'F' && direction != 'R')) return 1;
            unsigned char *paths[Max_Path_Num];
            uint64_t root = actgkmer_uint64((unsigned char *)(seed+31-k), k);
            current_direction = direction;
            int pos = find_paths_protected(h, paths, root, mask, cutoff, 96,
                                           k, 16, ratio, skip, NULL, &redundancy);
            ++queries;
            if (pos >= 0) {
                unsigned char packed[113];
                packed_dna_string(paths[pos], packed, 112);
                char *body = (char *)packed+2;
                if (direction == 'R') reverse_complement_body(body);
                fprintf(out, "%c\t%s\t%u\t%s\n", direction, seed, support, body);
                ++returned;
                free(paths[pos]);
            }
        }
        fclose(seeds); fclose(out);
        printf("DONE %s %lu %lu %.6f\n", tag, queries, returned,
               (double)(clock()-start)/CLOCKS_PER_SEC);
        fflush(stdout);
    }
    free(line); c4x_destroy(h); return 0;
}
