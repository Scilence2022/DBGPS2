#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include "ketopt.h"
#include "dbgps3_core.h"
#include "dbgps3_server.h"

int main(int argc, char *argv[])
{
    int c=1, k = 31, p = DBGPS3_KC_BITS, block_size = 10000000, n_thread = 3;
    int min_cov = 1, index_len = 16, payload_len = 140, path_len;
    double max_cov_ratio = 2.55;
    double max_r = 0;
    int mem_test = 0;
    int read_len = 200;
    int skip_ratios = 0;
    char *output_prefix = NULL;
    int server_mode = 0;
    int server_port = 9300;
    int max_clients = 10;

    unsigned int min_index = 101010102, max_index = 101295684, index_byte_num, max_cov = 0;
    uint64_t pr_uint64;

    char p1[] = "CCTGCAGAGTAGCATGTC";
    int pr_len = 18;

    ketopt_t o = KETOPT_INIT;
    while ((c = ketopt(&o, argc, argv, 1, "k:p:b:t:c:l:a:b:i:C:r:R:m:s:o:S::P:L:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'p') p = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 'c') min_cov = atoi(o.arg);
        else if (c == 'C') max_cov = atoi(o.arg);
        else if (c == 'a') min_index = atoi(o.arg);
        else if (c == 'b') max_index = atoi(o.arg);
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 'r') max_cov_ratio = atof(o.arg);
        else if (c == 'R') max_r = atof(o.arg);
        else if (c == 'm') mem_test = 1;
        else if (c == 's') skip_ratios = atoi(o.arg);
        else if (c == 'o') output_prefix = o.arg;
        else if (c == 'S') { server_mode = 1; if (o.arg) server_port = atoi(o.arg); }
        else if (c == 'P') server_port = atoi(o.arg);
        else if (c == 'L') max_clients = atoi(o.arg);
    }

    if (argc - o.ind < 1) {
        fprintf(stderr, "\n************************************************************************\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "**  DBGPS v2 - De Bruijn Graph based strand assembler for DNA data storage  **\n");
        fprintf(stderr, "**     Version 2.0 (2025-05-14)  Author: Lifu Song lifu.song@outlook.com      **\n");
        fprintf(stderr, "**                                                                    **\n");
        fprintf(stderr, "************************************************************************\n\n");
        fprintf(stderr, "Usage: DBGPS3 [options] <input file> \n");
        fprintf(stderr, "                       [Supporting formats: *fq, *fa, *fq.gz, *fa.gz]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size [%d]\n", k);
        fprintf(stderr, "  -i INT     length of index [%d] bp\n", index_len);
        fprintf(stderr, "  -l INT     data encoding length [%d] bp\n", payload_len);
        fprintf(stderr, "  -t INT     number of threads [%d]\n", n_thread);
        fprintf(stderr, "  -c INT     k-mer coverage cut-off for exclusion of noise k-mers [%d]\n", min_cov);
        fprintf(stderr, "  -C INT     Switch on k-mer coverage testing mode.  [%d]\n", max_cov);
        fprintf(stderr, "  -a INT     Initial index [%d]\n", min_index);
        fprintf(stderr, "  -b INT     End index [%d]\n", max_index);
        fprintf(stderr, "  -r FLOAT   Maximum coverage ratio [%.1f]\n", max_cov_ratio);
        fprintf(stderr, "  -R FLOAT   Upper bound for ratio range iteration (0 = disabled) [%.1f]\n", max_r);
        fprintf(stderr, "  -s INT     Number of initial positions to skip ratio checks for [%d]\n", skip_ratios);
        fprintf(stderr, "  -o STR     Output prefix for log file [none]\n");
        fprintf(stderr, "  -S [PORT]  Start in server mode (default port: %d)\n", server_port);
        fprintf(stderr, "  -P INT     Server TCP port [%d]\n", server_port);
        fprintf(stderr, "  -L INT     Server max concurrent clients [%d]\n", max_clients);
        fprintf(stderr, "\n");
        return 1;
    }
    if (p < DBGPS3_KC_BITS) {
        fprintf(stderr, "ERROR: -p should be at least %d\n", DBGPS3_KC_BITS);
        return 1;
    }

    if(max_cov < min_cov){max_cov = min_cov;}
    index_byte_num = index_len>>2;
    pr_uint64 = dbgps3_primer_uint64(p1, pr_len, index_len, k);
    path_len = payload_len;

    /* ---- Load the de Bruijn graph ---- */
    clock_t start, end;
    double kc_dur;

    fprintf(stderr, "Counting k-mers ......\n");
    start = clock();

    dbgps3_graph_t *g = dbgps3_load(k, p, n_thread, block_size, read_len,
                                      argc - o.ind, argv + o.ind);
    if (!g) {
        fprintf(stderr, "ERROR: failed to load input files\n");
        return 1;
    }

    /* Override graph parameters set during loading */
    g->index_len = index_len;
    g->index_byte_num = index_byte_num;
    g->path_len = path_len;

    end = clock();
    kc_dur = (double)(end - start);

    /* ---- Server mode ---- */
    if (server_mode) {
        fprintf(stderr, "\nK-mer counting time: %f seconds\n", (kc_dur/CLOCKS_PER_SEC));
        fprintf(stderr, "Starting DBGPS3 server on port %d (max clients: %d)...\n",
                server_port, max_clients);

        dbgps3_server_t srv;
        srv.graph = g;
        srv.port = server_port;
        srv.max_clients = max_clients;
        srv.running = 0;

        int rc = dbgps3_server_start(&srv);
        dbgps3_destroy(g);
        return rc;
    }

    /* ---- Batch CLI mode (original behavior) ---- */
    FILE *log_file = NULL;
    FILE *strands_file = NULL;
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

    #define PRINT_BOTH(format, ...) do { \
        fprintf(stdout, format, ##__VA_ARGS__); \
        if (log_file != NULL) fprintf(log_file, format, ##__VA_ARGS__); \
    } while(0)

    PRINT_BOTH("========== DBGPS Parameters ==========\n");
    PRINT_BOTH("k-mer size (-k): %d\n", k);
    PRINT_BOTH("Prefix length (-p): %d\n", p);
    PRINT_BOTH("Number of threads (-t): %d\n", n_thread);
    PRINT_BOTH("Minimum coverage cutoff (-c): %d\n", min_cov);
    PRINT_BOTH("Maximum coverage for testing (-C): %d\n", max_cov);
    PRINT_BOTH("Index length (-i): %d bp\n", index_len);
    PRINT_BOTH("Data encoding length (-l): %d bp\n", payload_len);
    PRINT_BOTH("Read length for k-mer counting: %d\n", read_len);
    PRINT_BOTH("Skip ratio checks for first (-s): %d positions\n", skip_ratios);
    PRINT_BOTH("Initial index (-a): %d\n", min_index);
    PRINT_BOTH("End index (-b): %d\n", max_index);
    PRINT_BOTH("Maximum coverage ratio (-r): %.2f\n", max_cov_ratio);
    PRINT_BOTH("Ratio range upper bound (-R): %.2f\n", max_r);
    PRINT_BOTH("Memory test mode: %s\n", mem_test ? "enabled" : "disabled");
    PRINT_BOTH("Path length (payload): %d\n", path_len);
    PRINT_BOTH("Index byte number: %d\n", index_byte_num);
    PRINT_BOTH("Output prefix: %s\n", output_prefix ? output_prefix : "none");
    PRINT_BOTH("======================================\n");

    PRINT_BOTH("\nK-mer counting time: %f seconds\n",(kc_dur/CLOCKS_PER_SEC));

    if(mem_test){
        fprintf(stderr, "Memory test mode activated! Press ENTER to continue...\n");
        getchar();
    }
    fprintf(stderr, "Decoding strands from index %d to ",min_index);
    fprintf(stderr, "%d ......\n",max_index);

    unsigned int id;
    unsigned int cov_cut;
    int c_index;
    int dec_strand_num = 0;
    double dec_dur;

    char result_buf[DBGPS3_MAX_PATH_LEN + 128]; /* enough for any single path result */

    if (max_r <= 0 || max_r <= max_cov_ratio) {
        /* Original behavior - use only the specified max_cov_ratio */
        for(cov_cut=min_cov; cov_cut<=max_cov; cov_cut++) {
            start = clock();
            dec_strand_num = 0;

            unsigned int eid;
            for(id=min_index; id<=max_index; id++) {
                fprintf(stderr, "Processing index: %u\t", id);
                eid = dbgps3_index_encode(id, index_len);
                c_index = dbgps3_find_path(g, pr_uint64 | (uint64_t)eid,
                                           cov_cut, max_cov_ratio, skip_ratios,
                                           result_buf, sizeof(result_buf));
                fprintf(stderr, "c_index: %u\n", (unsigned)c_index);

                if(c_index > -1) {
                    dec_strand_num++;
                    if (output_prefix != NULL && strands_file != NULL) {
                        fprintf(strands_file, "%d\t%s\n", id, result_buf);
                    } else {
                        printf("%d\t%s\n", id, result_buf);
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

                unsigned int eid;
                for(id=min_index; id<=max_index; id++) {
                    fprintf(stderr, "Processing index: %u\r", id);
                    eid = dbgps3_index_encode(id, index_len);
                    c_index = dbgps3_find_path(g, pr_uint64 | (uint64_t)eid,
                                               cov_cut, r, skip_ratios,
                                               result_buf, sizeof(result_buf));
                    if(c_index > -1) {
                        dec_strand_num++;
                        if (output_prefix != NULL && strands_file != NULL) {
                            fprintf(strands_file, "%d\t%s\n", id, result_buf);
                        } else {
                            printf("%d\t%s\n", id, result_buf);
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

    dbgps3_destroy(g);

    if (log_file != NULL) fclose(log_file);
    if (strands_file != NULL) fclose(strands_file);

    #undef PRINT_BOTH
    return 0;
}
