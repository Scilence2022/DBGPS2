/* dbgps2_profile: instrumented re-implementation of the DBGPS2 path search
 * used for the thermal-aging manuscript.
 *
 * It reuses the k-mer counter, coverage lookup, strand layout and CRC check of
 * DBGPS2.c verbatim (the file is #included), counts the reads once, and then:
 *   1. runs the bounded breadth-first path search for a grid of (c, r) values
 *      and classifies every index as decoded (correct / wrong), no path,
 *      path overflow, no CRC-valid path, or ambiguous;
 *   2. profiles the true-path coverage of every reference strand, including
 *      the coverage of the three alternative (noise) branches at every step;
 *   3. writes the k-mer coverage spectrum split into reference and noise k-mers.
 *
 * The search logic mirrors find_paths_protected() in DBGPS2.c; the decoded
 * counts are checked against the DBGPS2 binary in the manuscript pipeline.
 *
 * Build: cc -O3 -I.. -o dbgps2_profile dbgps2_profile.c ../dbgps_codec.c \
 *        ../dbgps_redundancy.c ../dbgps_strand.c ../kthread.c -lz -lpthread
 */
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>
#define main dbgps2_main
#include "DBGPS2.c"
#undef main

enum { OUT_OK = 0, OUT_WRONG, OUT_NOPATH, OUT_OVERFLOW, OUT_NOVALID, OUT_AMBIG, OUT_N };
static const char out_char[OUT_N] = { 'S', 'W', 'B', 'X', 'V', 'A' };

typedef struct {
    path_node *nodes;
    path_node **term, **next;
} search_buf_t;

typedef struct {
    kc_c4x_t *h;
    uint64_t mask, pr;
    int k, index_len, payload_len, path_len, seq_len, skip, cut;
    double r;
    const dbgps_codec_t *codec;
    const dbgps_redundancy_t *red;
    dbgps_strand_layout_t layout;
    uint64_t min_id, n_id;
    unsigned char **ref_packed; /* expected index+path packed DNA, NULL if absent */
    search_buf_t *bufs;
    unsigned char *outcome;   /* per id */
    uint32_t *nodes_used;     /* per id */
    uint16_t *width;          /* per id */
    unsigned char *accepted; /* truth-independent codewords, zero when rejected */
} grid_ctx_t;

static double wall(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static int search_one(grid_ctx_t *g, search_buf_t *b, uint64_t seed, unsigned char *out,
                      uint32_t *n_nodes, uint16_t *max_width)
{
    path_node *all = b->nodes;
    path_node **term = b->term, **next = b->next;
    int k = g->k, cut = g->cut;
    double r = g->r;
    all[0].pre_node = NULL;
    all[0].kmer = seed;
    all[0].cov = kmer_cov(min_hash_key(seed, k), g->mask, g->h);
    term[0] = &all[0];
    int term_num = 1, all_num = 1, p = 0, width = 1;
    path_node next_nodes[4];
    while (p <= g->path_len && term_num > 0) {
        int new_num = 0;
        for (int n_p = 0; n_p < term_num; ++n_p) {
            path_node *pn = term[n_p];
            int a_cov = pn->cov / r;
            double upper = pn->cov * r;
            int b_cov = upper > KC_MAX ? KC_MAX : (int)upper;
            if (a_cov < cut + 1) a_cov = cut + 1;
            for (int i = 0; i < 4; ++i) {
                next_nodes[i].kmer = ((pn->kmer << 2) + i) & g->mask;
                next_nodes[i].cov = kmer_cov(min_hash_key(next_nodes[i].kmer, k), g->mask, g->h);
                int eligible = p < g->skip ? (int)next_nodes[i].cov >= cut + 1
                    : (int)next_nodes[i].cov >= a_cov && (int)next_nodes[i].cov <= b_cov;
                if (eligible) {
                    if (new_num >= Max_Path_Num || all_num >= Max_Path_Num * Max_Path_Len) {
                        *n_nodes = all_num; *max_width = Max_Path_Num;
                        return OUT_OVERFLOW;
                    }
                    all[all_num].pre_node = pn;
                    all[all_num].kmer = next_nodes[i].kmer;
                    all[all_num].cov = next_nodes[i].cov;
                    next[new_num++] = &all[all_num++];
                }
            }
        }
        path_node **t = term; term = next; next = t;
        term_num = new_num;
        if (term_num > width) width = term_num;
        ++p;
    }
    *n_nodes = all_num; *max_width = width > 65535 ? 65535 : width;
    if (term_num == 0) return OUT_NOPATH;
    unsigned char candidate[(Max_Path_Len + 3) / 4];
    unsigned char corrected[sizeof(candidate)], raw[sizeof(candidate)];
    int found = 0;
    for (int i = 0; i < term_num; ++i) {
        if (g->path_len % 4 == 0 && g->index_len % 4 == 0)
            nodes_path(term[i], candidate, g->path_len, k, g->index_len / 4);
        else nodes_path_packed(term[i], candidate, g->path_len, g->index_len);
        if (dbgps_strand_check(g->codec, g->red, &g->layout, candidate, corrected, raw) < 0)
            continue;
        if (found) {
            if (!memcmp(out, corrected, g->seq_len)) continue;
            return OUT_AMBIG;
        }
        memcpy(out, corrected, g->seq_len);
        found = 1;
    }
    return found ? OUT_OK : OUT_NOVALID;
}

static void grid_worker(void *data, long i, int tid)
{
    grid_ctx_t *g = (grid_ctx_t*)data;
    uint64_t id = g->min_id + (uint64_t)i, eid;
    if (g->codec) dbgps_codec_encode_index(g->codec, id, g->index_len, &eid);
    else eid = index_encode(id, g->index_len);
    unsigned char out[(Max_Path_Len + 3) / 4];
    uint32_t nn = 0; uint16_t w = 0;
    int o = search_one(g, &g->bufs[tid], g->pr | eid, out, &nn, &w);
    memset(g->accepted + (size_t)i*g->seq_len, 0, g->seq_len);
    if (o == OUT_OK) memcpy(g->accepted + (size_t)i*g->seq_len, out, g->seq_len);
    if (o == OUT_OK && g->ref_packed[i] && memcmp(out, g->ref_packed[i], g->seq_len)) o = OUT_WRONG;
    g->outcome[i] = (unsigned char)o;
    g->nodes_used[i] = nn;
    g->width[i] = w;
}

/* pack ACGT text into canonical 2-bit MSB-first bytes */
static void pack_dna(const char *s, int n, unsigned char *out)
{
    memset(out, 0, (n + 3) / 4);
    for (int i = 0; i < n; ++i)
        out[i / 4] |= (seq_nt4_table[(uint8_t)s[i]] & 3) << (6 - 2 * (i % 4));
}

static uint64_t kmer_of(const char *s, int k)
{
    uint64_t x = 0;
    for (int i = 0; i < k; ++i) x = x << 2 | (seq_nt4_table[(uint8_t)s[i]] & 3);
    return x;
}

static int parse_list(const char *s, double *v, int max)
{
    int n = 0;
    while (*s && n < max) {
        char *e;
        v[n++] = strtod(s, &e);
        if (e == s) return -1;
        s = *e == ',' ? e + 1 : e;
    }
    return n;
}

static long peak_rss_mb(void)
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return ru.ru_maxrss / (1024 * 1024);
#else
    return ru.ru_maxrss / 1024;
#endif
}

#define RH_BINS 400 /* log2 ratio bins from -10 to +10, width 0.05 */
static int rbin(double lr)
{
    int b = (int)floor((lr + 10.0) / 0.05);
    return b < 0 ? 0 : b >= RH_BINS ? RH_BINS - 1 : b;
}

int main(int argc, char *argv[])
{
    int k = 31, index_len = 16, payload_len = 140, n_thread = 8, skip = 0, F = 0;
    const char *ref_fn = NULL, *prefix = "profile", *map_fn = NULL, *red_spec = "crc16";
    const char *c_list = "1,2,3,4,5", *r_list = "0,2.55";
    uint64_t min_id = 101010102, max_id = 101295684;
    int do_stats = 1, validate_every = 0;
    const char *trace_fn = NULL; /* ids whose per-step coverages are dumped */
    ketopt_t o = KETOPT_INIT;
    int c;
    while ((c = ketopt(&o, argc, argv, 1, "k:i:l:t:s:x:o:M:E:c:r:a:b:F:NT:v:", 0)) >= 0) {
        if (c == 'k') k = atoi(o.arg);
        else if (c == 'i') index_len = atoi(o.arg);
        else if (c == 'l') payload_len = atoi(o.arg);
        else if (c == 't') n_thread = atoi(o.arg);
        else if (c == 's') skip = atoi(o.arg);
        else if (c == 'x') ref_fn = o.arg;
        else if (c == 'o') prefix = o.arg;
        else if (c == 'M') map_fn = o.arg;
        else if (c == 'E') red_spec = o.arg;
        else if (c == 'c') c_list = o.arg;
        else if (c == 'r') r_list = o.arg;
        else if (c == 'a') min_id = strtoull(o.arg, NULL, 10);
        else if (c == 'b') max_id = strtoull(o.arg, NULL, 10);
        else if (c == 'F') F = atoi(o.arg);
        else if (c == 'N') do_stats = 0;
        else if (c == 'T') trace_fn = o.arg;
        else if (c == 'v') validate_every = atoi(o.arg);
    }
    if (argc - o.ind < 1 || !ref_fn) {
        fprintf(stderr, "Usage: dbgps2_profile -x ref.fa [-M map] [-c 1,2,..] [-r 0,2.55,..] [-o prefix] reads...\n"
                        "       r = 0 disables the ratio constraint (original DBGPS behaviour)\n");
        return 1;
    }
    dbgps_codec_t mapping = {0};
    const dbgps_codec_t *codec = NULL;
    char err[512];
    if (map_fn) {
        if (dbgps_codec_load(&mapping, map_fn, err, sizeof err) < 0) { fprintf(stderr, "%s\n", err); return 1; }
        codec = &mapping;
    }
    dbgps_redundancy_t red = {0};
    if (dbgps_redundancy_init(&red, red_spec, err, sizeof err) < 0) { fprintf(stderr, "%s\n", err); return 1; }
    grid_ctx_t g = {0};
    if (dbgps_strand_layout(codec, &red, index_len, payload_len, &g.layout) < 0) return 1;
    g.k = k; g.index_len = index_len; g.payload_len = payload_len; g.skip = skip;
    g.path_len = g.layout.path_bases;
    g.seq_len = (g.layout.total_bases + 3) / 4;
    g.mask = (1ULL << k * 2) - 1;
    char p1[] = "CCTGCAGAGTAGCATGTC";
    int pr_len = 18;
    g.pr = primer_uint64(p1, pr_len, index_len, k);
    g.codec = codec; g.red = &red;
    g.min_id = min_id; g.n_id = max_id - min_id + 1;
    int read_len = 2 * pr_len + index_len + g.path_len;
    if (read_len < 200) read_len = 200;

    /* ---- reference strands ---- */
    g.ref_packed = (unsigned char**)calloc(g.n_id, sizeof(*g.ref_packed));
    char **ref_seq = (char**)calloc(g.n_id, sizeof(*ref_seq));
    {
        gzFile fp = gzopen(ref_fn, "r");
        if (!fp) { fprintf(stderr, "cannot open %s\n", ref_fn); return 1; }
        kseq_t *ks = kseq_init(fp);
        int n = 0;
        while (kseq_read(ks) >= 0) {
            uint64_t id = strtoull(ks->name.s, NULL, 10);
            if (id < min_id || id > max_id) continue;
            if ((int)ks->seq.l < pr_len + index_len + g.path_len) continue;
            uint64_t i = id - min_id;
            ref_seq[i] = strdup(ks->seq.s);
            g.ref_packed[i] = (unsigned char*)malloc(g.seq_len);
            pack_dna(ks->seq.s + pr_len, index_len + g.path_len, g.ref_packed[i]);
            ++n;
        }
        kseq_destroy(ks); gzclose(fp);
        fprintf(stderr, "Loaded %d reference strands\n", n);
    }

    /* ---- count reads ---- */
    double t0 = wall();
    kc_c4x_t *h;
    if (F) h = count_files_filtered(argc - o.ind, argv + o.ind, k, KC_BITS, 10000000, n_thread, read_len, F);
    else {
        h = count_file(argv[o.ind], k, KC_BITS, 10000000, n_thread, read_len);
        for (int f = o.ind + 1; h && f < argc; ++f)
            if (!count_file2(argv[f], h, k, KC_BITS, 10000000, n_thread, read_len)) h = NULL;
    }
    if (!h) { fprintf(stderr, "counting failed\n"); return 1; }
    double t_count = wall() - t0;
    uint64_t n_distinct = 0;
    for (int i = 0; i < 1 << KC_BITS; ++i) n_distinct += kh_size(h->h[i]);
    g.h = h;
    char fn[2048];
    FILE *fo;
    snprintf(fn, sizeof fn, "%s.info.tsv", prefix);
    fo = fopen(fn, "w");
    fprintf(fo, "count_seconds\t%.2f\ndistinct_kmers\t%" PRIu64 "\npeak_rss_mb_after_count\t%ld\n",
            t_count, n_distinct, peak_rss_mb());
    fclose(fo);
    fprintf(stderr, "Counted %" PRIu64 " distinct k-mers in %.1f s\n", n_distinct, t_count);


    /* Oracle diagnostics are never used for seed enumeration or path decisions. */
    const int nstep = g.path_len + 1;
    uint16_t *truth = calloc((size_t)g.n_id * (nstep + 1), sizeof(uint16_t));
    uint16_t *alts = calloc((size_t)g.n_id * nstep * 3, sizeof(uint16_t));
    assert(truth && alts);
    double mean_cov_all = 0, cv_all = 0, mincov_all = 0;
    uint64_t nref_profile = 0;
    int start_pos = pr_len + index_len - k;
    for (uint64_t i = 0; i < g.n_id; ++i) {
        if (!ref_seq[i]) continue;
        uint16_t *v = truth + (size_t)i * (nstep + 1);
        uint64_t parent = kmer_of(ref_seq[i] + start_pos, k);
        v[0] = kmer_cov(min_hash_key(parent,k),g.mask,h);
        double su=0, sq=0; int minv=KC_MAX;
        for (int p=0;p<nstep;++p) {
            uint64_t child=kmer_of(ref_seq[i]+start_pos+p+1,k);
            v[p+1]=kmer_cov(min_hash_key(child,k),g.mask,h);
            if(p<nstep-1){su+=v[p+1]; sq+=(double)v[p+1]*v[p+1];
                if(v[p+1]<minv)minv=v[p+1];}
            int j=0;
            for(int a=0;a<4;++a){
                uint64_t alt=((parent<<2)|a)&g.mask;
                if(alt!=child)alts[((size_t)i*nstep+p)*3+j++]=kmer_cov(min_hash_key(alt,k),g.mask,h);
            }
            parent=child;
        }
        double m=su/(nstep-1), var=sq/(nstep-1)-m*m;
        mean_cov_all+=m; cv_all+=m>0?sqrt(var>0?var:0)/m:0; mincov_all+=minv;
        ++nref_profile;
    }
    /* Compact, lossless per-step evidence, little-endian uint16 on this host. */
    snprintf(fn,sizeof fn,"%s.truth.u16",prefix);fo=fopen(fn,"wb");
    fwrite(truth,sizeof(uint16_t),(size_t)g.n_id*(nstep+1),fo);fclose(fo);
    snprintf(fn,sizeof fn,"%s.alts.u16",prefix);fo=fopen(fn,"wb");
    fwrite(alts,sizeof(uint16_t),(size_t)g.n_id*nstep*3,fo);fclose(fo);
    snprintf(fn,sizeof fn,"%s.oracle.json",prefix);fo=fopen(fn,"w");
    fprintf(fo,"{\"min_id\":%" PRIu64 ",\"n_id\":%" PRIu64 ",\"steps\":%d,\"dtype\":\"little-endian uint16\",\"k\":%d}\n",min_id,g.n_id,nstep,k);fclose(fo);

    /* ---- coverage spectrum: reference vs noise k-mers ---- */
    if (do_stats) {
        kc_c4x_t *rt = count_file(ref_fn, k, KC_BITS, 10000000, n_thread, 400);
        uint64_t *spec_ref = (uint64_t*)calloc(KC_MAX + 1, 8), *spec_noise = (uint64_t*)calloc(KC_MAX + 1, 8);
        uint64_t ref_total = 0, ref_absent = 0;
        for (int i = 0; i < 1 << KC_BITS; ++i) {
            kc_c4_t *t = h->h[i], *r = rt->h[i];
            for (khint_t j = 0; j < kh_end(t); ++j) {
                if (!kh_exist(t, j)) continue;
                uint64_t key = kh_key(t, j);
                int cnt = key & KC_MAX;
                khint_t x = kc_c4_get(r, key >> KC_BITS << KC_BITS);
                if (x != kh_end(r)) spec_ref[cnt]++; else spec_noise[cnt]++;
            }
            for (khint_t j = 0; j < kh_end(r); ++j) {
                if (!kh_exist(r, j)) continue;
                ref_total++;
                if (kc_c4_get(t, kh_key(r, j) >> KC_BITS << KC_BITS) == kh_end(t)) ref_absent++;
            }
        }
        c4x_destroy(rt);
        snprintf(fn, sizeof fn, "%s.spectrum.tsv", prefix);
        fo = fopen(fn, "w");
        fprintf(fo, "# ref_kmers_total=%" PRIu64 " ref_kmers_absent=%" PRIu64 "\ncov\tref\tnoise\n", ref_total, ref_absent);
        for (int i = 0; i <= KC_MAX; ++i)
            if (spec_ref[i] || spec_noise[i]) fprintf(fo, "%d\t%" PRIu64 "\t%" PRIu64 "\n", i, spec_ref[i], spec_noise[i]);
        fclose(fo);
        free(spec_ref); free(spec_noise);

        /* ---- true-path profile ---- */
        uint64_t rh_true[RH_BINS] = {0}, rh_alt1[RH_BINS] = {0}, rh_alt2[RH_BINS] = {0};
        int steps = g.path_len + 1;
        double *pos_sum = (double*)calloc(steps, sizeof(double)), *pos_alt = (double*)calloc(steps, sizeof(double));
        uint64_t *pos_zero = (uint64_t*)calloc(steps, 8);
        uint64_t n_prof = 0;
        snprintf(fn, sizeof fn, "%s.strandstats.tsv", prefix);
        fo = fopen(fn, "w");
        fprintf(fo, "id\tseed_cov\tmin_cov\tmedian_cov\tmean_cov\tr_needed\tr_needed_skip5\tmax_alt_cov\tmax_alt_frac\tn_alt_ge2\tn_alt_ge2_r255\tcv\n");
        int start = pr_len + index_len - k;
        for (uint64_t i = 0; i < g.n_id; ++i) {
            const char *s = ref_seq[i];
            if (!s) continue;
            n_prof++;
            uint64_t parent = kmer_of(s + start, k);
            int pc = kmer_cov(min_hash_key(parent, k), g.mask, h), seed = pc;
            int minc = KC_MAX, maxalt = 0, nalt2 = 0, nalt2r = 0;
            double rneed = 1.0, rneed5 = 1.0, sum = 0, sum2 = 0, maxaltfrac = 0;
            int covs[Max_Path_Len];
            for (int p = 0; p < steps; ++p) {
                uint64_t child = kmer_of(s + start + 1 + p, k);
                int tc = kmer_cov(min_hash_key(child, k), g.mask, h);
                covs[p] = tc;
                sum += tc; sum2 += (double)tc * tc;
                pos_sum[p] += tc;
                if (tc == 0) pos_zero[p]++;
                if (tc < minc) minc = tc;
                double rr = (pc > 0 && tc > 0) ? (pc > tc ? (double)pc / tc : (double)tc / pc) : 1e9;
                if (rr > rneed) rneed = rr;
                if (p >= 5 && rr > rneed5) rneed5 = rr;
                if (pc > 0 && tc > 0) rh_true[rbin(log2((double)tc / pc))]++;
                for (int a = 0; a < 4; ++a) {
                    uint64_t alt = ((parent << 2) | a) & g.mask;
                    if (alt == child) continue;
                    int ac = kmer_cov(min_hash_key(alt, k), g.mask, h);
                    if (ac <= 0) continue;
                    pos_alt[p] += 1;
                    if (pc > 0) {
                        rh_alt1[rbin(log2((double)ac / pc))]++;
                        if (ac >= 2) rh_alt2[rbin(log2((double)ac / pc))]++;
                    }
                    if (ac > maxalt) maxalt = ac;
                    if (pc > 0 && (double)ac / pc > maxaltfrac) maxaltfrac = (double)ac / pc;
                    if (ac >= 2) {
                        nalt2++;
                        if (pc > 0 && ac >= (int)(pc / 2.55) && ac <= pc * 2.55) nalt2r++;
                    }
                }
                parent = child; pc = tc;
            }
            /* median */
            int tmp[Max_Path_Len];
            memcpy(tmp, covs, steps * sizeof(int));
            for (int a = 1; a < steps; ++a) { int v = tmp[a], b = a - 1; while (b >= 0 && tmp[b] > v) { tmp[b + 1] = tmp[b]; --b; } tmp[b + 1] = v; }
            double mean = sum / steps, var = sum2 / steps - mean * mean;
            fprintf(fo, "%" PRIu64 "\t%d\t%d\t%d\t%.2f\t%.4g\t%.4g\t%d\t%.4g\t%d\t%d\t%.4f\n", min_id + i, seed, minc,
                    tmp[steps / 2], mean, rneed, rneed5, maxalt, maxaltfrac, nalt2, nalt2r,
                    mean > 0 ? sqrt(var > 0 ? var : 0) / mean : 0);
        }
        fclose(fo);
        snprintf(fn, sizeof fn, "%s.ratiohist.tsv", prefix);
        fo = fopen(fn, "w");
        fprintf(fo, "log2_ratio_lo\ttrue\talt_ge1\talt_ge2\n");
        for (int b = 0; b < RH_BINS; ++b)
            fprintf(fo, "%.2f\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n", -10.0 + b * 0.05, rh_true[b], rh_alt1[b], rh_alt2[b]);
        fclose(fo);
        snprintf(fn, sizeof fn, "%s.posprofile.tsv", prefix);
        fo = fopen(fn, "w");
        fprintf(fo, "step\tmean_cov\tfrac_zero\tmean_alt_branches\n");
        for (int p = 0; p < steps; ++p)
            fprintf(fo, "%d\t%.3f\t%.6f\t%.4f\n", p, pos_sum[p] / n_prof, (double)pos_zero[p] / n_prof, pos_alt[p] / n_prof);
        fclose(fo);
        free(pos_sum); free(pos_alt); free(pos_zero);
        fprintf(stderr, "Profiled %" PRIu64 " strands\n", n_prof);
    }

    /* ---- optional per-step traces: true child and the three alternatives ---- */
    if (trace_fn) {
        FILE *ft = fopen(trace_fn, "r");
        snprintf(fn, sizeof fn, "%s.trace.tsv", prefix);
        fo = fopen(fn, "w");
        fprintf(fo, "id\tstep\tparent_cov\ttrue_cov\talt1\talt2\talt3\n");
        uint64_t id;
        int start = pr_len + index_len - k;
        while (ft && fscanf(ft, "%" SCNu64, &id) == 1) {
            if (id < min_id || id > max_id || !ref_seq[id - min_id]) continue;
            const char *s = ref_seq[id - min_id];
            uint64_t parent = kmer_of(s + start, k);
            int pc = kmer_cov(min_hash_key(parent, k), g.mask, h);
            for (int p = 0; p <= g.path_len; ++p) {
                uint64_t child = kmer_of(s + start + 1 + p, k);
                int tc = kmer_cov(min_hash_key(child, k), g.mask, h), alt[3], na = 0;
                for (int a = 0; a < 4; ++a) {
                    uint64_t x = ((parent << 2) | a) & g.mask;
                    if (x != child) alt[na++] = kmer_cov(min_hash_key(x, k), g.mask, h);
                }
                fprintf(fo, "%" PRIu64 "\t%d\t%d\t%d\t%d\t%d\t%d\n", id, p, pc, tc, alt[0], alt[1], alt[2]);
                parent = child; pc = tc;
            }
        }
        if (ft) fclose(ft);
        fclose(fo);
    }

    /* ---- decoding grid ---- */
    double cs[256], rs[256];
    int nc = parse_list(c_list, cs, 256), nr = parse_list(r_list, rs, 256);
    if (nc < 1 || nr < 1) { fprintf(stderr, "bad -c/-r list\n"); return 1; }
    g.bufs = (search_buf_t*)calloc(n_thread, sizeof(search_buf_t));
    for (int t = 0; t < n_thread; ++t) {
        g.bufs[t].nodes = (path_node*)malloc(sizeof(path_node) * Max_Path_Num * Max_Path_Len);
        g.bufs[t].term = (path_node**)malloc(sizeof(path_node*) * Max_Path_Num);
        g.bufs[t].next = (path_node**)malloc(sizeof(path_node*) * Max_Path_Num);
    }
    g.accepted = (unsigned char*)calloc((size_t)g.n_id, g.seq_len);
    g.outcome = (unsigned char*)malloc(g.n_id);
    g.nodes_used = (uint32_t*)malloc(4 * g.n_id);
    g.width = (uint16_t*)malloc(2 * g.n_id);
    int ncombo = nc * nr;
    unsigned char *mat = (unsigned char*)malloc((size_t)ncombo * g.n_id);
    snprintf(fn, sizeof fn, "%s.grid.tsv", prefix);
    FILE *fg = fopen(fn, "w");
    fprintf(fg, "c\tr\tskip\tn_ref\tdecoded_ref\tcorrect\twrong\tno_path\toverflow\tno_valid\tambiguous\t"
                "spurious_nonref\tmean_nodes_ref\tmean_width_ref\tseconds\tceil_cov\tceil_end\tceil_seed\tceil_adm\tceil_cont\tlookahead_losses\tseed_losses\tratio_losses\tsearch_gap\toverflow_adm\tambiguous_adm\twrong_adm\tno_valid_adm\tnoise_per_admitted_step\tmean_true_cov\tmean_cov_cv\tmean_mincov\tnative_checks\n");
    snprintf(fn,sizeof fn,"%s.accepted.u8",prefix); FILE *faccepted=fopen(fn,"wb");
    snprintf(fn,sizeof fn,"%s.blind_status.u8",prefix); FILE *fblind=fopen(fn,"wb");
    assert(faccepted && fblind);
    int combo = 0;
    for (int ri = 0; ri < nr; ++ri) for (int ci = 0; ci < nc; ++ci, ++combo) {
        g.cut = (int)cs[ci];
        g.r = rs[ri] <= 0 ? 1e6 : rs[ri];
        double ts = wall();
        kt_for(n_thread, grid_worker, &g, (long)g.n_id);
        double el = wall() - ts;
        assert(fwrite(g.accepted, g.seq_len, g.n_id, faccepted)==g.n_id);
        for(uint64_t bi=0;bi<g.n_id;++bi) fputc(g.outcome[bi]==OUT_WRONG?OUT_OK:g.outcome[bi],fblind);
        uint64_t cnt[OUT_N] = {0}, nref = 0, spur = 0;
        double nodes = 0, width = 0;
        for (uint64_t i = 0; i < g.n_id; ++i) {
            if (g.ref_packed[i]) {
                nref++; cnt[g.outcome[i]]++;
                nodes += g.nodes_used[i]; width += g.width[i];
            } else if (g.outcome[i] == OUT_OK) spur++;
            mat[(size_t)combo * g.n_id + i] = out_char[g.outcome[i]];
        }

        uint64_t uc=0,ue=0,us=0,ua=0,cont=0,oa=0,aa=0,wa=0,va=0,noise=0,nchecked=0;
        for(uint64_t i=0;i<g.n_id;++i){
            if(!g.ref_packed[i])continue;
            uint16_t *v=truth+(size_t)i*(nstep+1);
            int covok=1,endok=0,adm=1,con=1; uint64_t nb=0;
            for(int p=0;p<nstep;++p){
                int a=v[p], b=v[p+1];
                int lo=(int)(a/g.r), hi=a*g.r>KC_MAX?KC_MAX:(int)(a*g.r);
                if(lo<g.cut+1)lo=g.cut+1;
                if(p<nstep-1){
                    if(b<g.cut+1)covok=0;
                    if(b<g.cut+1 || (p>=g.skip && (b<lo || b>hi)))adm=0;
                    if(b<g.cut+1 || (p>=g.skip && (!a || !b || (double)b>a*g.r || (double)b<a/g.r)))con=0;
                    for(int j=0;j<3;++j){int ac=alts[((size_t)i*nstep+p)*3+j];
                        if(p<g.skip ? ac>=g.cut+1 : ac>=lo && ac<=hi)++nb;
                    }
                }else{
                    /* Any terminal successor completes the same reported CRC codeword. */
                    int last_adm=0,last_cont=0;
                    for(int j=-1;j<3;++j){int ac=j<0?b:alts[((size_t)i*nstep+p)*3+j];
                        if(ac>=g.cut+1){
                            endok=1;
                            if(p<g.skip || (ac>=lo && ac<=hi))last_adm=1;
                            if(p<g.skip || (a && (double)ac<=a*g.r && (double)ac>=a/g.r))last_cont=1;
                        }
                    }
                    if(!last_adm)adm=0;if(!last_cont)con=0;
                }
            }
            endok=endok&&covok;int seedok=endok&&(v[0]>0 || g.skip>0);
            uc+=covok;ue+=endok;us+=seedok;ua+=adm;cont+=con;
            if(adm){noise+=nb;oa+=g.outcome[i]==OUT_OVERFLOW;aa+=g.outcome[i]==OUT_AMBIG;wa+=g.outcome[i]==OUT_WRONG;va+=g.outcome[i]==OUT_NOVALID;}
            assert(!adm || seedok);assert(!con || adm);
            assert(g.outcome[i]!=OUT_OK || adm);
            if(validate_every>0 && i%(uint64_t)validate_every==0){
                uint64_t eid; if(codec)dbgps_codec_encode_index(codec,min_id+i,index_len,&eid);else eid=index_encode(min_id+i,index_len);
                unsigned char out[(Max_Path_Len+3)/4];uint32_t nn;uint16_t width;
                int inst=search_one(&g,&g.bufs[0],g.pr|eid,out,&nn,&width);
                unsigned char *paths[Max_Path_Num]={0};
                int native=find_paths_protected(h,paths,g.pr|eid,g.mask,g.cut,payload_len,k,index_len,g.r,skip,codec,&red);
                assert((native>=0)==(inst==OUT_OK));
                if(native>=0){assert(!memcmp(paths[native],out,g.seq_len));free(paths[native]);}
                ++nchecked;
            }
        }
        assert(cnt[OUT_OK]<=ua && cont<=ua && ua<=us && us<=ue && ue<=uc);
        fprintf(fg, "%d\t%g\t%d\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
                    "\t%" PRIu64 "\t%" PRIu64 "\t%.1f\t%.2f\t%.2f",
                g.cut, rs[ri], skip, nref, cnt[OUT_OK] + cnt[OUT_WRONG], cnt[OUT_OK], cnt[OUT_WRONG], cnt[OUT_NOPATH],
                cnt[OUT_OVERFLOW], cnt[OUT_NOVALID], cnt[OUT_AMBIG], spur, nodes / nref, width / nref, el);

        fprintf(fg,"\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.9g\t%.9g\t%.9g\t%.9g\t%" PRIu64 "\n",
            uc,ue,us,ua,cont,uc-ue,ue-us,us-ua,ua-cnt[OUT_OK],oa,aa,wa,va,ua?(double)noise/ua/(nstep-1):0,
            mean_cov_all/nref_profile,cv_all/nref_profile,mincov_all/nref_profile,nchecked);
        fflush(fg);
        fprintf(stderr, "c=%d r=%g decoded %" PRIu64 " (%.1fs)\n", g.cut, rs[ri], cnt[OUT_OK], el);
    }
    fclose(fg); fclose(faccepted); fclose(fblind);
    snprintf(fn,sizeof fn,"%s.blind_meta.json",prefix);fo=fopen(fn,"w");
    fprintf(fo,"{\"min_id\":%" PRIu64 ",\"n_id\":%" PRIu64 ",\"seq_bytes\":%d,\"n_combos\":%d,\"ordering\":\"r outer, c inner; lists supplied in command\",\"accepted_status\":0}\n",min_id,g.n_id,g.seq_len,ncombo);fclose(fo);
    snprintf(fn, sizeof fn, "%s.outcomes.tsv", prefix);
    fo = fopen(fn, "w");
    fprintf(fo, "id");
    for (int ri = 0; ri < nr; ++ri) for (int ci = 0; ci < nc; ++ci) fprintf(fo, "\tc%d_r%g", (int)cs[ci], rs[ri]);
    fprintf(fo, "\n");
    for (uint64_t i = 0; i < g.n_id; ++i) {
        if (!g.ref_packed[i]) continue;
        fprintf(fo, "%" PRIu64 "\t", min_id + i);
        for (int cb = 0; cb < ncombo; ++cb) fputc(mat[(size_t)cb * g.n_id + i], fo);
        fputc('\n', fo);
    }
    fclose(fo);
    snprintf(fn, sizeof fn, "%s.info.tsv", prefix);
    fo = fopen(fn, "a");
    fprintf(fo, "peak_rss_mb_final\t%ld\n", peak_rss_mb());
    fclose(fo);
    return 0;
}
