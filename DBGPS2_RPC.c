#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include "ketopt.h"
#include "dbgps_bloom.h"
#include "dbgps_core.h"
#include "dbgps_rpc_server.h"

/* ketopt implements required short-option arguments only. Consume -S's
   optional port explicitly so -S never swallows the input filename or -P. */
static int optional_server_port(ketopt_t *o, int argc, char **argv, int *port)
{
    uint64_t value;
    int attached = o->pos != 0;
    if (o->i >= argc) return 0;
    const char *text = argv[o->i] + (attached ? o->pos : 0);
    if (dbgps_parse_index(text, &value) < 0) return attached ? -1 : 0;
    if (value < 1 || value > 65535) return -1;
    *port = (int)value;
    if (o->n_args) ketopt_permute(argv, o->i, o->n_args);
    ++o->i;
    o->pos = 0;
    o->ind = o->i - o->n_args;
    return 0;
}

int main(int argc, char *argv[])
{
    int c=1, k = 31, p = DBGPS_KC_BITS, block_size = 10000000, n_thread = 3;
    int min_cov = 1, index_len = 16, payload_len = 140, path_len;
    double max_cov_ratio = 2.55;
    double max_r = 0;
    int mem_test = 0;
    int bloom_bits_per_key = 0;
    int count_bloom_shift = 0;
    int read_len = 200;
    int skip_ratios = 0;
    char *output_prefix = NULL;
    int server_mode = 0;
    int server_port = 9300;
    int max_clients = 10;

    uint64_t min_index = 101010102, max_index = 101295684;
    int index_byte_num;
    unsigned int max_cov = 0;
    const char *mapping_file = NULL;
    const char *redundancy_spec = DBGPS_DEFAULT_REDUNDANCY;
    uint64_t pr_uint64;

    char p1[] = "CCTGCAGAGTAGCATGTC";
    int pr_len = 18;

    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:c:l:a:b:i:C:r:R:m:s:o:M:E:B:F:SP:L:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'C') max_cov = atoi(o.arg);
        else if (c == 'a' || c == 'b') {
            if (dbgps_parse_index(o.arg, c == 'a' ? &min_index : &max_index) < 0) {
                fprintf(stderr, "ERROR: -%c requires an unsigned decimal index\n", c);
                return 1;
            }
        }
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 'r') max_cov_ratio = atof(o.arg);
        else if (c == 'R') max_r = atof(o.arg);
        else if (c == 'm') mem_test = 1;
        else if (c == 's') skip_ratios = atoi(o.arg);
        else if (c == 'o') output_prefix = o.arg;
        else if (c == 'M') mapping_file = o.arg;
        else if (c == 'E') redundancy_spec = o.arg;
        else if (c == 'F') {
            if (dbgps_bloom_parse_shift(o.arg, &count_bloom_shift) < 0) {
                fprintf(stderr, "ERROR: -F requires 0 (disabled) or log2 bitmap bits from 23 to 40\n");
                return 1;
            }
        }
        else if (c == 'B') {
            if (dbgps_bloom_parse_bits(o.arg, &bloom_bits_per_key) < 0) {
                fprintf(stderr, "ERROR: -B requires an integer from 0 to 64 (0 disables Bloom filtering)\n");
                return 1;
            }
        }
        else if (c == 'S') {
            server_mode = 1;
            if (optional_server_port(&o, argc, argv, &server_port) < 0) {
                fprintf(stderr, "ERROR: server port must be an integer from 1 to 65535\n");
                return 1;
            }
        }
        else if (c == 'P') server_port = atoi(o.arg);
        else if (c == 'L') max_clients = atoi(o.arg);
        else {
            fprintf(stderr, "ERROR: unknown option or missing option value\n");
            return 1;
        }
    }

    if (argc - o.ind < 1) {
        fprintf(stderr, "\n************************************************************************\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "**  DBGPS2-RPC - De Bruijn Graph based strand assembler for DNA data storage  **\n");
        fprintf(stderr, "**     Version 2.0 (2025-05-14)  Author: Lifu Song lifu.song@outlook.com      **\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "************************************************************************\n\n");
        fprintf(stderr, "Usage: DBGPS2-RPC [options] <input file> \n");
        fprintf(stderr, "                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
        fprintf(stderr, "  -i INT     length of index [%d] bp\n", index_len);
        fprintf(stderr, "  -l INT     data encoding length [%d] bp\n", payload_len);
        fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
        fprintf(stderr, "  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [%d]\n", min_cov);
        fprintf(stderr, "  -C INT     Switch on k-mer coverage testing mode.  [%d]\n", max_cov);
        fprintf(stderr, "  -a INT     Initial index [%"PRIu64"]\n", min_index);
        fprintf(stderr, "  -b INT     End index [%"PRIu64"]\n", max_index);
        fprintf(stderr, "  -r FLOAT   Maximum coverage ratio [%.1f]\n", max_cov_ratio);
        fprintf(stderr, "  -R FLOAT   Upper bound for ratio range iteration (0 = disabled) [%.1f]\n", max_r);
        fprintf(stderr, "  -s INT     Number of initial positions to skip ratio checks for [%d]\n", skip_ratios);
        fprintf(stderr, "  -o STR     Output prefix for log file [none]\n");
        fprintf(stderr, "  -M FILE    Custom DNA/bit mapping [none; built-in A=00,C=01,G=10,T=11]\n");
        fprintf(stderr, "  -E CODE    Redundancy: " DBGPS_REDUNDANCY_OPTIONS " [" DBGPS_DEFAULT_REDUNDANCY "]\n");
        fprintf(stderr, "  -B INT     Bloom query prefilter bits/key (1..64; 0 disables) [0]\n");
        fprintf(stderr, "  -F INT     Two-pass singleton filter log2 bits (23..40; 0 disables) [0]\n");
        fprintf(stderr, "  -S [PORT]  Start in server mode (default port: %d)\n", server_port);
        fprintf(stderr, "  -P INT     Server TCP port [%d]\n", server_port);
        fprintf(stderr, "  -L INT     Server max concurrent clients [%d]\n", max_clients);
        fprintf(stderr, "\n");
        return 1;
    }
    if (p != DBGPS_KC_BITS) {
        fprintf(stderr, "ERROR: -p must be %d\n", DBGPS_KC_BITS);
        return 1;
    }

    dbgps_graph_t config = {0};
    config.k = k;
    config.index_len = index_len;
    config.path_len = payload_len;
    char mapping_error[512];
    if (dbgps_set_mapping(&config, mapping_file, mapping_error, sizeof(mapping_error)) < 0) {
        fprintf(stderr, "ERROR: %s\n", mapping_error);
        return 1;
    }
    if (dbgps_set_redundancy(&config, redundancy_spec, mapping_error, sizeof(mapping_error)) < 0) {
        fprintf(stderr, "ERROR: %s\n", mapping_error);
        dbgps_codec_destroy(&config.codec);
        return 1;
    }
    char redundancy_label[DBGPS_REDUNDANCY_SPEC_SIZE];
    dbgps_redundancy_format(&config.redundancy, redundancy_label, sizeof(redundancy_label));
    dbgps_strand_layout_t layout;
    dbgps_strand_layout(config.codec.bases ? &config.codec : NULL, &config.redundancy,
                         index_len, payload_len, &layout);
    uint64_t encoded_max;
    if (dbgps_validate_layout(&config, index_len, payload_len, &index_byte_num) < 0 ||
        k - index_len > pr_len ||
        (!server_mode && (min_index > max_index ||
                         dbgps_encode_index(&config, max_index, index_len, &encoded_max) < 0))) {
        fprintf(stderr, "ERROR: incompatible strand layout or index range; index and payload must\n"
                        "       contain whole mapping blocks and whole bytes, index must fit k,\n"
                        "       k-index must fit the primer, total strand length must be < %d bp,\n"
                        "       and redundancy must fit its codeword capacity (RS: data + parity <= 255 bytes)\n", DBGPS_MAX_PATH_LEN);
        dbgps_codec_destroy(&config.codec);
        return 1;
    }
    if(max_cov < min_cov){max_cov = min_cov;}
    pr_uint64 = dbgps_primer_uint64(p1, pr_len, index_len, k);
    path_len = payload_len;
    if (config.codec.bases || 2 * pr_len + (int)layout.total_bases > read_len)
        read_len = 2 * pr_len + layout.total_bases;

    /* ---- Batch output, opened and reported before counting as in DBGPS2 ---- */
    FILE *log_file = NULL;
    FILE *strands_file = NULL;
    int exit_status = 0;

    #define PRINT_BOTH(format, ...) do { \
        fprintf(stdout, format, ##__VA_ARGS__); \
        if (log_file != NULL) fprintf(log_file, format, ##__VA_ARGS__); \
    } while(0)

    if (!server_mode) {
        if (output_prefix != NULL) {
            char log_filename[1024];
            sprintf(log_filename, "%s.log", output_prefix);
            log_file = fopen(log_filename, "w");
            if (log_file == NULL) {
                fprintf(stderr, "Warning: Could not open log file %s for writing\n", log_filename);
            } else {
                fprintf(stderr, "Writing parameters to log file: %s\n", log_filename);
            }

            char strands_filename[1024];
            sprintf(strands_filename, "%s.strands", output_prefix);
            strands_file = fopen(strands_filename, "w");
            if (strands_file == NULL) {
                fprintf(stderr, "Warning: Could not open strands file %s for writing\n", strands_filename);
            } else {
                fprintf(stderr, "Will write assembled strands to: %s\n", strands_filename);
            }
        }

        PRINT_BOTH("========== DBGPS Parameters ==========\n");
        PRINT_BOTH("k-mer size (-k): %d\n", k);
        PRINT_BOTH("Prefix length (-p): %d\n", p);
        PRINT_BOTH("Number of threads (-t): %d\n", n_thread);
        PRINT_BOTH("Bloom counting filter (-F): %d (0 disables; singletons removed when enabled)\n", count_bloom_shift);
        PRINT_BOTH("Minimum coverage cutoff (-c): %d\n", min_cov);
        PRINT_BOTH("Maximum coverage for testing (-C): %d\n", max_cov);
        PRINT_BOTH("Index length (-i): %d bp\n", index_len);
        PRINT_BOTH("Data encoding length (-l): %d bp\n", payload_len);
        PRINT_BOTH("Read length for k-mer counting: %d\n", read_len);
        PRINT_BOTH("Skip ratio checks for first (-s): %d positions\n", skip_ratios);
        PRINT_BOTH("Initial index (-a): %"PRIu64"\n", min_index);
        PRINT_BOTH("End index (-b): %"PRIu64"\n", max_index);
        PRINT_BOTH("Maximum coverage ratio (-r): %.2f\n", max_cov_ratio);
        PRINT_BOTH("Ratio range upper bound (-R): %.2f\n", max_r);
        PRINT_BOTH("Memory test mode: %s\n", mem_test ? "enabled" : "disabled");
        PRINT_BOTH("Path length (payload + redundancy): %u\n", layout.path_bases);
        PRINT_BOTH("Index byte number: %d\n", index_byte_num);
        if (mapping_file) {
            PRINT_BOTH("Mapping: %s (%u bases <-> %u bits; %s)\n", mapping_file,
                       config.codec.bases, config.codec.bits,
                       config.codec.byte_fast ? "byte lookup" : "block lookup");
        } else {
            PRINT_BOTH("Mapping: none (built-in A=00,C=01,G=10,T=11; direct 2-bit packing)\n");
        }
        PRINT_BOTH("Redundancy (-E): %s; %u parity bytes; %u bp; zero padding: %u bits\n",
                   redundancy_label, layout.parity_bytes, layout.parity_bases, layout.padding_bits);
        PRINT_BOTH("Output prefix: %s\n", output_prefix ? output_prefix : "none");
        PRINT_BOTH("======================================\n");
    }

    /* ---- Load the de Bruijn graph ---- */
    clock_t start, end;
    double kc_dur;

    fprintf(stderr, "Counting k-mers ......\n");
    start = clock();

    dbgps_graph_t *g = dbgps_load_filtered(k, p, n_thread, block_size, read_len,
                                      argc - o.ind, argv + o.ind, count_bloom_shift);
    if (!g) {
        fprintf(stderr, "ERROR: failed to load input files or allocate graph/filter memory\n");
        dbgps_codec_destroy(&config.codec);
        exit_status = 1;
        goto cleanup;
    }
    g->codec = config.codec; /* Transfer ownership to the graph. */
    g->redundancy = config.redundancy;
    if (dbgps_set_bloom(g, bloom_bits_per_key) < 0) {
        fprintf(stderr, "ERROR: failed to allocate Bloom query filter\n");
        exit_status = 1;
        goto cleanup;
    }

    /* Override graph parameters set during loading */
    g->index_len = index_len;
    g->index_byte_num = index_byte_num;
    g->path_len = path_len;

    if (!server_mode)
        PRINT_BOTH("Bloom query filter (-B): %d bits/key; %zu bitmap bytes\n",
                   bloom_bits_per_key, dbgps_bloom_bytes(g));

    end = clock();
    kc_dur = (double)(end - start);

    /* ---- Server mode ---- */
    if (server_mode) {
        if (bloom_bits_per_key)
            fprintf(stderr, "Bloom query filter: %d bits/key; %zu bitmap bytes\n",
                    bloom_bits_per_key, dbgps_bloom_bytes(g));
        if (count_bloom_shift)
            fprintf(stderr, "Bloom counting filter: 2^%d bits; singletons removed\n", count_bloom_shift);
        if (mapping_file)
            fprintf(stderr, "Mapping: %s (%u bases <-> %u bits)\n", mapping_file, g->codec.bases, g->codec.bits);
        else fprintf(stderr, "Mapping: none (built-in A=00,C=01,G=10,T=11; direct 2-bit packing)\n");
        fprintf(stderr, "Redundancy (-E): %s; %u parity bytes; %u bp; zero padding: %u bits\n",
                redundancy_label, layout.parity_bytes, layout.parity_bases, layout.padding_bits);
        fprintf(stderr, "\nK-mer counting time: %f seconds\n", (kc_dur/CLOCKS_PER_SEC));
        fprintf(stderr, "Starting DBGPS2-RPC server on port %d (max clients: %d)...\n",
                server_port, max_clients);

        dbgps_rpc_server_t srv;
        srv.graph = g;
        srv.port = server_port;
        srv.max_clients = max_clients;
        srv.running = 0;

        int rc = dbgps_rpc_server_start(&srv);
        dbgps_destroy(g);
        return rc;
    }

    /* ---- Batch CLI mode ---- */
    PRINT_BOTH("\nK-mer counting time: %f seconds\n",(kc_dur/CLOCKS_PER_SEC));

    if(mem_test){
        fprintf(stderr, "Memory test mode activated! Press ENTER to continue...\n");
        getchar();
    }
    fprintf(stderr, "Decoding strands from index %"PRIu64" to ",min_index);
    fprintf(stderr, "%"PRIu64" ......\n",max_index);

    uint64_t id;
    unsigned int cov_cut;
    int c_index;
    int dec_strand_num = 0;
    double dec_dur;

    char result_buf[DBGPS_MAX_PATH_LEN + 128]; /* enough for any single path result */

    if (max_r <= 0 || max_r <= max_cov_ratio) {
        /* Original behavior - use only the specified max_cov_ratio */
        for(cov_cut=min_cov; cov_cut<=max_cov; cov_cut++) {
            start = clock();
            dec_strand_num = 0;

            uint64_t eid;
            for(id=min_index; id<=max_index; id++) {
                if (dbgps_encode_index(g, id, index_len, &eid) < 0) continue;
                c_index = dbgps_find_path(g, pr_uint64 | (uint64_t)eid,
                                           cov_cut, max_cov_ratio, skip_ratios,
                                           path_len, index_byte_num,
                                           result_buf, sizeof(result_buf));

                if(c_index > -1) {
                    dec_strand_num++;
                    if (output_prefix != NULL && strands_file != NULL) {
                        fprintf(strands_file, "%"PRIu64"\t%s\n", id, result_buf);
                    } else {
                        printf("%"PRIu64"\t%s\n", id, result_buf);
                    }
                }
            }

            end = clock();
            dec_dur = (double)(end - start);
            PRINT_BOTH("\nCoverage cutoff: %d, Ratio: %.2f", cov_cut, max_cov_ratio);
            PRINT_BOTH("\nStrands decoded: %d ", dec_strand_num);
            PRINT_BOTH("\nStrand decoding time: %f seconds\n",(dec_dur/CLOCKS_PER_SEC));
        }
    } else {
        /* New behavior - iterate through different r values */
        for(double r = max_cov_ratio; r <= max_r; r += 0.05) {
            for(cov_cut=min_cov; cov_cut<=max_cov; cov_cut++) {
                start = clock();
                dec_strand_num = 0;

                uint64_t eid;
                for(id=min_index; id<=max_index; id++) {
                    if (dbgps_encode_index(g, id, index_len, &eid) < 0) continue;
                    c_index = dbgps_find_path(g, pr_uint64 | (uint64_t)eid,
                                               cov_cut, r, skip_ratios,
                                               path_len, index_byte_num,
                                               result_buf, sizeof(result_buf));
                    if(c_index > -1) {
                        dec_strand_num++;
                        if (output_prefix != NULL && strands_file != NULL) {
                            fprintf(strands_file, "%"PRIu64"\t%s\n", id, result_buf);
                        } else {
                            printf("%"PRIu64"\t%s\n", id, result_buf);
                        }
                    }
                }

                end = clock();
                dec_dur = (double)(end - start);
                PRINT_BOTH("\nCoverage cutoff: %d, Ratio: %.2f", cov_cut, r);
                PRINT_BOTH("\nStrands decoded: %d ", dec_strand_num);
                PRINT_BOTH("\nStrand decoding time: %f seconds\n",(dec_dur/CLOCKS_PER_SEC));
            }
        }
    }

cleanup:
    dbgps_destroy(g);

    if (log_file != NULL) fclose(log_file);
    if (strands_file != NULL) fclose(strands_file);

    #undef PRINT_BOTH
    return exit_status;
}
