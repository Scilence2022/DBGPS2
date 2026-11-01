/* Thin experiment adapter: the graph and path search are the unmodified DBGPS2
 * implementation at the repository revision recorded in provenance.json.
 * Input seeds come from sequencing reads, never from reference sequences.
 * Build from repository root; see run_experiment.py. */
#define main dbgps2_original_main
#include "../../../DBGPS2.c"
#undef main

int main(int argc, char **argv)
{
    if (argc != 5 && argc != 6) {
        fprintf(stderr, "Usage: %s reads.fa[.gz] seeds.tsv output_prefix threads [include_no_ratio]\n", argv[0]);
        return 1;
    }
    const int k = 31, index_len = 16, payload_len = 96;
    const uint64_t mask = (1ULL << 62) - 1;
    clock_t started = clock();
    kc_c4x_t *h = count_file(argv[1], k, KC_BITS, 10000000, atoi(argv[4]), 200);
    if (!h) { fprintf(stderr, "Cannot count input\n"); return 1; }
    fprintf(stderr, "graph_cpu_seconds %.6f\n", (double)(clock() - started) / CLOCKS_PER_SEC);
    dbgps_redundancy_t redundancy = {0};
    if (dbgps_redundancy_init(&redundancy, "none", NULL, 0) < 0) return 1;
    const char *configs[] = {"c0_r2.55", "c1_r2.55", "c0_no_ratio"};
    const int cutoffs[] = {0, 1, 0}, skips[] = {0, 0, 97};
    int config_count = argc == 6 && atoi(argv[5]) ? 3 : 2;
    for (int cfg = 0; cfg < config_count; ++cfg) {
        FILE *seeds = fopen(argv[2], "r");
        char file[4096];
        snprintf(file, sizeof(file), "%s.%s.tsv", argv[3], configs[cfg]);
        FILE *out = fopen(file, "w");
        if (!seeds || !out) return 1;
        fprintf(out, "seed\tsupport\tpayload\n");
        char seed[64]; unsigned support;
        unsigned long attempted = 0, returned = 0;
        started = clock();
        while (fscanf(seeds, "%63s %u", seed, &support) == 2) {
            if (strlen(seed) != 31 || !support) return 1;
            unsigned char *paths[Max_Path_Num];
            uint64_t start = actgkmer_uint64((unsigned char *)seed, k);
            int pos = find_paths_protected(h, paths, start, mask, cutoffs[cfg], payload_len,
                k, index_len, 2.55, skips[cfg], NULL, &redundancy);
            ++attempted;
            if (pos >= 0) {
                unsigned char dna[113];
                packed_dna_string(paths[pos], dna, 112);
                /* Output is two real primer bases + the complete 110-base body.
                 * No synthetic sequence, padding, or reference bases were added. */
                fprintf(out, "%s\t%u\t%s\n", seed, support, dna + 2);
                ++returned;
                free(paths[pos]);
            }
        }
        fprintf(stderr, "%s queries %lu returned %lu search_cpu_seconds %.6f\n",
            configs[cfg], attempted, returned, (double)(clock() - started) / CLOCKS_PER_SEC);
        fclose(seeds); fclose(out);
    }
    c4x_destroy(h);
    return 0;
}
