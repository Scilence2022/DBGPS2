#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <math.h>
#include <limits.h>
#include <inttypes.h>
#include "dbgps3_core.h"
#include "dbgps3_server.h"
#include "minjson.h"

/* ---- JSON-RPC 2.0 Error Codes ---- */
#define JRPC_PARSE_ERROR      -32700
#define JRPC_INVALID_REQUEST  -32600
#define JRPC_METHOD_NOT_FOUND -32601
#define JRPC_INVALID_PARAMS   -32602
#define JRPC_INTERNAL_ERROR   -32603

/* ---- Helper: build a JSON-RPC error response ---- */
static mj_val jrpc_error(int code, const char *message, mj_val id)
{
    mj_val resp = mj_object();
    mj_object_set(&resp, "jsonrpc", mj_string("2.0"));
    mj_object_set(&resp, "id", id);

    mj_val err = mj_object();
    mj_object_set(&err, "code", mj_number(code));
    mj_object_set(&err, "message", mj_string(message));

    mj_object_set(&resp, "error", err);
    return resp;
}

/* ---- Helper: build a JSON-RPC success response ---- */
static mj_val jrpc_result(mj_val result, mj_val id)
{
    mj_val resp = mj_object();
    mj_object_set(&resp, "jsonrpc", mj_string("2.0"));
    mj_object_set(&resp, "id", id);
    mj_object_set(&resp, "result", result);
    return resp;
}

/* ---- Handler: kmer/cov ---- */
static mj_val handle_kmer_cov(dbgps3_graph_t *g, const mj_val *params, mj_val id)
{
    mj_val pv;
    pv = mj_object_get(params, "kmer");
    const char *kmer_str = mj_get_string(&pv, NULL);
    if (!kmer_str) {
        return jrpc_error(JRPC_INVALID_PARAMS, "Missing 'kmer' parameter", id);
    }

    int k = dbgps3_get_k(g);
    int len = (int)strlen(kmer_str);
    if (len != k) {
        char msg[128];
        snprintf(msg, sizeof(msg), "kmer length %d != k=%d", len, k);
        return jrpc_error(JRPC_INVALID_PARAMS, msg, id);
    }

    uint64_t kmer_key = dbgps3_actgkmer_hashkey(kmer_str, k);
    unsigned short cov = dbgps3_kmer_cov(g, kmer_key);

    mj_val result = mj_object();
    mj_object_set(&result, "coverage", mj_number(cov));
    mj_object_set(&result, "kmer", mj_string(kmer_str));
    return jrpc_result(result, id);
}

/* ---- Handler: path/find ---- */
static int optional_int(const mj_val *params, const char *key, int fallback,
                         int minimum, int maximum, int *value)
{
    mj_val field = mj_object_get(params, key);
    if (!mj_object_has(params, key)) {
        *value = fallback;
        return 0;
    }
    if (field.type != MJ_NUMBER || !isfinite(field.number) ||
        field.number < minimum || field.number > maximum ||
        (double)(int)field.number != field.number) return -1;
    *value = (int)field.number;
    return 0;
}

typedef struct {
    int cut_cov, skip_ratios, index_len, index_bytes, path_len;
    double ratio;
} path_options_t;

static int read_path_options(const dbgps3_graph_t *g, const mj_val *params,
                              int raw_kmer, path_options_t *out)
{
    if (optional_int(params, "cutCov", 1, 0, DBGPS3_KC_MAX - 1, &out->cut_cov) < 0 ||
        optional_int(params, "skipRatios", 0, 0, INT_MAX, &out->skip_ratios) < 0 ||
        optional_int(params, "pathLen", g->path_len, 0, DBGPS3_MAX_PATH_LEN, &out->path_len) < 0)
        return -1;
    out->index_len = 0;
    if (!raw_kmer && optional_int(params, "indexLen", g->index_len, 0, g->k, &out->index_len) < 0)
        return -1;
    mj_val ratio = mj_object_get(params, "maxCovRatio");
    if (mj_object_has(params, "maxCovRatio") && ratio.type != MJ_NUMBER) return -1;
    out->ratio = mj_get_number(&ratio, 2.55);
    if (!isfinite(out->ratio) || out->ratio < 1 ||
        (!raw_kmer && g->k - out->index_len > 18)) return -1;
    return dbgps3_validate_layout(g, out->index_len, out->path_len, &out->index_bytes);
}

static mj_val handle_path_find(dbgps3_graph_t *g, const mj_val *params, mj_val id)
{
    mj_val pv = mj_object_get(params, "index");
    uint64_t index;
    if (pv.type == MJ_STRING) {
        if (dbgps3_parse_index(pv.string, &index) < 0)
            return jrpc_error(JRPC_INVALID_PARAMS, "Invalid unsigned decimal 'index'", id);
    } else if (pv.type == MJ_NUMBER && isfinite(pv.number) && pv.number >= 0 &&
               pv.number <= 9007199254740991.0 && (double)(uint64_t)pv.number == pv.number) {
        index = (uint64_t)pv.number;
    } else {
        return jrpc_error(JRPC_INVALID_PARAMS, "Missing or invalid 'index' parameter", id);
    }

    path_options_t options;
    uint64_t eid;
    if (read_path_options(g, params, 0, &options) < 0 ||
        dbgps3_encode_index(g, index, options.index_len, &eid) < 0)
        return jrpc_error(JRPC_INVALID_PARAMS, "Incompatible mapping layout or index capacity", id);

    /* Default primer: same as CLI */
    const char *primer = "CCTGCAGAGTAGCATGTC";
    int pr_len = 18;
    uint64_t pr_uint64 = dbgps3_primer_uint64(primer, pr_len, options.index_len, g->k);
    uint64_t start_kmer = pr_uint64 | eid;

    char result_buf[DBGPS3_MAX_PATH_LEN + 128];
    int rc = dbgps3_find_path(g, start_kmer, options.cut_cov, options.ratio, options.skip_ratios,
                               options.path_len, options.index_bytes,
                               result_buf, sizeof(result_buf));

    mj_val result = mj_object();
    if (pv.type == MJ_STRING) {
        char decimal[32];
        snprintf(decimal, sizeof(decimal), "%"PRIu64, index);
        mj_object_set(&result, "index", mj_string(decimal));
    } else mj_object_set(&result, "index", mj_number((double)index));
    if (rc == 0) {
        mj_object_set(&result, "found", mj_bool(1));
        mj_object_set(&result, "sequence", mj_string(result_buf));
    } else {
        mj_object_set(&result, "found", mj_bool(0));
    }
    return jrpc_result(result, id);
}

/* ---- Handler: path/findKmer ---- */
static mj_val handle_path_find_kmer(dbgps3_graph_t *g, const mj_val *params, mj_val id)
{
    mj_val pv;
    pv = mj_object_get(params, "kmer");
    const char *kmer_str = mj_get_string(&pv, NULL);
    if (!kmer_str) {
        return jrpc_error(JRPC_INVALID_PARAMS, "Missing 'kmer' parameter", id);
    }

    int k = dbgps3_get_k(g);
    int len = (int)strlen(kmer_str);
    if (len != k) {
        char msg[128];
        snprintf(msg, sizeof(msg), "kmer length %d != k=%d", len, k);
        return jrpc_error(JRPC_INVALID_PARAMS, msg, id);
    }

    path_options_t options;
    if (read_path_options(g, params, 1, &options) < 0)
        return jrpc_error(JRPC_INVALID_PARAMS, "Incompatible mapping layout or path parameters", id);

    uint64_t start_kmer = dbgps3_actgkmer_uint64(kmer_str, k);

    char result_buf[DBGPS3_MAX_PATH_LEN + 128];
    int rc = dbgps3_find_path(g, start_kmer, options.cut_cov, options.ratio, options.skip_ratios,
                               options.path_len, 0,
                               result_buf, sizeof(result_buf));

    mj_val result = mj_object();
    mj_object_set(&result, "kmer", mj_string(kmer_str));
    if (rc == 0) {
        mj_object_set(&result, "found", mj_bool(1));
        mj_object_set(&result, "sequence", mj_string(result_buf));
    } else {
        mj_object_set(&result, "found", mj_bool(0));
    }
    return jrpc_result(result, id);
}

/* ---- Handler: graph/info ---- */
static mj_val handle_graph_info(dbgps3_graph_t *g, const mj_val *params, mj_val id)
{
    (void)params;
    mj_val result = mj_object();
    mj_object_set(&result, "k", mj_number(dbgps3_get_k(g)));
    mj_object_set(&result, "index_len", mj_number(dbgps3_get_index_len(g)));
    mj_object_set(&result, "index_byte_num", mj_number(dbgps3_get_index_byte_num(g)));
    mj_object_set(&result, "path_len", mj_number(dbgps3_get_path_len(g)));
    mj_object_set(&result, "kmer_count", mj_number(dbgps3_get_kmer_count(g)));
    mj_val mapping = mj_object();
    int custom = g->codec.bases != 0;
    mj_object_set(&mapping, "custom", mj_bool(custom));
    mj_object_set(&mapping, "bases", mj_number(custom ? g->codec.bases : 1));
    mj_object_set(&mapping, "bits", mj_number(custom ? g->codec.bits : 2));
    mj_object_set(&mapping, "mode", mj_string(!custom ? "canonical" :
                  g->codec.byte_fast ? "byte_lookup" : "block_lookup"));
    mj_object_set(&result, "mapping", mapping);
    return jrpc_result(result, id);
}

/* ---- Handler: graph/load ---- */
static mj_val handle_graph_load(dbgps3_graph_t *g, const mj_val *params, mj_val id)
{
    mj_val files_val = mj_object_get(params, "files");
    if (files_val.type != MJ_ARRAY || files_val.array.count < 1) {
        return jrpc_error(JRPC_INVALID_PARAMS, "Missing or invalid 'files' parameter", id);
    }

    int n = files_val.array.count;
    char **file_paths = (char**)malloc(n * sizeof(char*));
    int i;
    for (i = 0; i < n; i++) {
        mj_val fv = mj_array_get(&files_val, i);
        file_paths[i] = (char*)mj_get_string(&fv, "");
    }

    int rc = dbgps3_load_more(g, n, file_paths);
    free(file_paths);

    if (rc != 0) {
        return jrpc_error(JRPC_INTERNAL_ERROR, "Failed to load files", id);
    }

    mj_val result = mj_object();
    mj_object_set(&result, "status", mj_string("ok"));
    return jrpc_result(result, id);
}

/* ---- Dispatch table ---- */
typedef struct {
    const char *method;
    mj_val (*handler)(dbgps3_graph_t*, const mj_val*, mj_val);
} method_entry;

static method_entry method_table[] = {
    {"kmer/cov",        handle_kmer_cov},
    {"path/find",       handle_path_find},
    {"path/findKmer",   handle_path_find_kmer},
    {"graph/info",      handle_graph_info},
    {"graph/load",      handle_graph_load},
    {NULL, NULL}
};

/* ---- Process a single JSON-RPC request ---- */
static mj_val process_request(dbgps3_server_t *srv, const mj_val *req)
{
    /* Validate request structure */
    if (req->type != MJ_OBJECT) {
        return jrpc_error(JRPC_INVALID_REQUEST, "Request must be a JSON object", mj_null());
    }

    /* Check jsonrpc version */
    mj_val version_val = mj_object_get(req, "jsonrpc");
    const char *version = mj_get_string(&version_val, NULL);
    if (!version || strcmp(version, "2.0") != 0) {
        return jrpc_error(JRPC_INVALID_REQUEST, "Invalid or missing jsonrpc version", mj_null());
    }

    /* Get method */
    mj_val method_val = mj_object_get(req, "method");
    const char *method = mj_get_string(&method_val, NULL);
    if (!method) {
        return jrpc_error(JRPC_INVALID_REQUEST, "Missing method", mj_null());
    }

    /* Get id (may be null) */
    mj_val id = mj_object_get(req, "id");

    /* Get params */
    mj_val params = mj_object_get(req, "params");
    if (params.type != MJ_OBJECT) {
        /* Default to empty object */
        params = mj_object();
    }

    /* Check for server/shutdown */
    if (strcmp(method, "server/shutdown") == 0) {
        mj_val result = mj_object();
        mj_object_set(&result, "status", mj_string("shutting down"));
        mj_val resp = jrpc_result(result, id);
        dbgps3_server_stop(srv);
        return resp;
    }

    /* Dispatch to handler */
    int i;
    for (i = 0; method_table[i].method != NULL; i++) {
        if (strcmp(method, method_table[i].method) == 0) {
            return method_table[i].handler(srv->graph, &params, id);
        }
    }

    return jrpc_error(JRPC_METHOD_NOT_FOUND, "Method not found", id);
}

/* ---- Line reader for socket ---- */
static int read_line(int fd, char *buf, int buf_size)
{
    int i = 0;
    while (i < buf_size - 1) {
        char c;
        ssize_t n = read(fd, &c, 1);
        if (n <= 0) return -1; /* connection closed or error */
        if (c == '\n') break;
        if (c != '\r') buf[i++] = c;
    }
    buf[i] = '\0';
    return i;
}

/* ---- Write all bytes to socket ---- */
static int write_all(int fd, const char *buf, int len)
{
    int written = 0;
    while (written < len) {
        ssize_t n = write(fd, buf + written, len - written);
        if (n <= 0) return -1;
        written += (int)n;
    }
    return written;
}

/* ---- Client handler thread ---- */
typedef struct {
    int fd;
    struct sockaddr_in addr;
    dbgps3_server_t *srv;
} client_info_t;

static void *client_handler(void *arg)
{
    client_info_t *info = (client_info_t*)arg;
    int fd = info->fd;
    dbgps3_server_t *srv = info->srv;

    char line[65536]; /* large buffer for JSON-RPC requests */

    while (srv->running) {
        int len = read_line(fd, line, sizeof(line));
        if (len < 0) break; /* connection closed */
        if (len == 0) continue; /* empty line */

        /* Parse JSON-RPC request */
        mj_val req = mj_parse(line);

        /* Process request */
        mj_val resp = process_request(srv, &req);

        /* Serialize and send response */
        char *resp_str = mj_serialize(&resp);
        if (resp_str) {
            int resp_len = (int)strlen(resp_str);
            write_all(fd, resp_str, resp_len);
            write_all(fd, "\n", 1);
            free(resp_str);
        }

        mj_free(&req);
        mj_free(&resp);

        /* Check if we just handled a shutdown */
        if (!srv->running) break;
    }

    close(fd);
    free(info);
    return NULL;
}

/* ---- Global server pointer for signal handler ---- */
static dbgps3_server_t *g_srv = NULL;

static void sigint_handler(int sig)
{
    (void)sig;
    if (g_srv) {
        g_srv->running = 0;
    }
}

/* ---- Server start ---- */
int dbgps3_server_start(dbgps3_server_t *srv)
{
    if (!srv || !srv->graph) return -1;

    g_srv = srv;
    srv->running = 1;

    /* Install SIGINT handler for graceful shutdown */
    struct sigaction sa;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = sigint_handler;
    sigaction(SIGINT, &sa, NULL);

    /* Create TCP socket */
    int server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd < 0) {
        perror("socket");
        return -1;
    }

    /* Allow address reuse */
    int opt = 1;
    setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

    /* Bind */
    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(srv->port);

    if (bind(server_fd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("bind");
        close(server_fd);
        return -1;
    }

    /* Listen */
    if (listen(server_fd, srv->max_clients) < 0) {
        perror("listen");
        close(server_fd);
        return -1;
    }

    fprintf(stderr, "[dbgps3_server] Listening on port %d (max clients: %d)\n",
            srv->port, srv->max_clients);
    fprintf(stderr, "[dbgps3_server] Graph loaded: k=%d, kmer_count=%ld\n",
            dbgps3_get_k(srv->graph), dbgps3_get_kmer_count(srv->graph));
    fprintf(stderr, "[dbgps3_server] Press Ctrl+C or send 'server/shutdown' to stop\n");

    /* Accept loop */
    while (srv->running) {
        struct sockaddr_in client_addr;
        socklen_t client_len = sizeof(client_addr);
        int client_fd = accept(server_fd, (struct sockaddr*)&client_addr, &client_len);

        if (client_fd < 0) {
            if (errno == EINTR) continue;
            perror("accept");
            break;
        }

        /* Create client info */
        client_info_t *info = (client_info_t*)malloc(sizeof(client_info_t));
        if (!info) {
            close(client_fd);
            continue;
        }
        info->fd = client_fd;
        info->addr = client_addr;
        info->srv = srv;

        /* Spawn client handler thread */
        pthread_t tid;
        if (pthread_create(&tid, NULL, client_handler, info) != 0) {
            perror("pthread_create");
            close(client_fd);
            free(info);
            continue;
        }
        pthread_detach(tid);

        char client_ip[INET_ADDRSTRLEN];
        inet_ntop(AF_INET, &client_addr.sin_addr, client_ip, sizeof(client_ip));
        fprintf(stderr, "[dbgps3_server] Client connected: %s:%d\n",
                client_ip, ntohs(client_addr.sin_port));
    }

    fprintf(stderr, "[dbgps3_server] Shutting down...\n");
    close(server_fd);
    g_srv = NULL;
    return 0;
}

void dbgps3_server_stop(dbgps3_server_t *srv)
{
    if (srv) {
        srv->running = 0;
        fprintf(stderr, "[dbgps3_server] Stop signal received\n");
    }
}
