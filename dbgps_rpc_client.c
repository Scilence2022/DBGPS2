/*
 * dbgps2-rpc-client - Interactive CLI client for DBGPS2-RPC server
 *
 * Connects to a DBGPS2-RPC server via TCP and sends JSON-RPC 2.0 requests.
 * Supports both interactive (REPL) and one-shot command modes.
 *
 * Usage:
 *   dbgps2-rpc-client [-h host] [-p port]                   # interactive mode
 *   echo '{"jsonrpc":"2.0","method":"graph/info",...}' | dbgps2-rpc-client -h host -p port
 *   dbgps2-rpc-client -h host -p port -c 'graph/info'      # shortcut command mode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "minjson.h"

#define DEFAULT_HOST "127.0.0.1"
#define DEFAULT_PORT 9300
#define BUF_SIZE 65536

static int connect_to_server(const char *host, int port)
{
    int fd = socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) {
        perror("socket");
        return -1;
    }

    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);

    if (inet_pton(AF_INET, host, &addr.sin_addr) <= 0) {
        fprintf(stderr, "Invalid address: %s\n", host);
        close(fd);
        return -1;
    }

    if (connect(fd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("connect");
        close(fd);
        return -1;
    }

    return fd;
}

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

static int read_line(int fd, char *buf, int buf_size)
{
    int i = 0;
    while (i < buf_size - 1) {
        char c;
        ssize_t n = read(fd, &c, 1);
        if (n <= 0) return -1;
        if (c == '\n') break;
        if (c != '\r') buf[i++] = c;
    }
    buf[i] = '\0';
    return i;
}

/* Build a JSON-RPC 2.0 request from method and params */
static char *build_request(const char *method, mj_val params, int id)
{
    mj_val req = mj_object();
    mj_object_set(&req, "jsonrpc", mj_string("2.0"));
    mj_object_set(&req, "method", mj_string(method));
    mj_object_set(&req, "id", mj_number(id));
    if (params.type != MJ_NULL) {
        mj_object_set(&req, "params", params);
    }

    char *str = mj_serialize(&req);
    mj_free(&req);
    return str;
}

/* Send a request and print the response */
static int send_request(int fd, const char *json_str, int pretty)
{
    /* Append newline */
    int len = (int)strlen(json_str);
    write_all(fd, json_str, len);
    write_all(fd, "\n", 1);

    /* Read response */
    char resp_buf[BUF_SIZE];
    int resp_len = read_line(fd, resp_buf, sizeof(resp_buf));
    if (resp_len < 0) {
        fprintf(stderr, "Error: connection closed\n");
        return -1;
    }

    if (pretty) {
        mj_val resp = mj_parse(resp_buf);
        char *pretty_str = mj_serialize_pretty(&resp);
        if (pretty_str) {
            printf("%s\n", pretty_str);
            free(pretty_str);
        }
        mj_free(&resp);
    } else {
        printf("%s\n", resp_buf);
    }

    return 0;
}

/* Shortcut: build and send a command by name */
static int send_shortcut(int fd, const char *cmd, int id)
{
    mj_val params = mj_object();
    char *request_str = NULL;

    if (strcmp(cmd, "info") == 0 || strcmp(cmd, "graph/info") == 0) {
        request_str = build_request("graph/info", params, id);
    }
    else if (strncmp(cmd, "cov ", 4) == 0) {
        mj_object_set(&params, "kmer", mj_string(cmd + 4));
        request_str = build_request("kmer/cov", params, id);
    }
    else if (strncmp(cmd, "find ", 5) == 0) {
        mj_object_set(&params, "index", mj_string(cmd + 5));
        request_str = build_request("path/find", params, id);
    }
    else if (strncmp(cmd, "findkmer ", 9) == 0) {
        mj_object_set(&params, "kmer", mj_string(cmd + 9));
        request_str = build_request("path/findKmer", params, id);
    }
    else if (strcmp(cmd, "shutdown") == 0 || strcmp(cmd, "server/shutdown") == 0) {
        request_str = build_request("server/shutdown", params, id);
    }
    else {
        fprintf(stderr, "Unknown shortcut: %s\n", cmd);
        fprintf(stderr, "Shortcuts: info, cov <kmer>, find <index>, findkmer <kmer>, shutdown\n");
        mj_free(&params);
        return -1;
    }

    int rc = send_request(fd, request_str, 1);
    free(request_str);
    return rc;
}

static void print_usage(void)
{
    fprintf(stderr, "dbgps2-rpc-client - Interactive CLI for DBGPS2-RPC server\n\n");
    fprintf(stderr, "Usage: dbgps2-rpc-client [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h HOST   Server hostname (default: %s)\n", DEFAULT_HOST);
    fprintf(stderr, "  -p PORT   Server port (default: %d)\n", DEFAULT_PORT);
    fprintf(stderr, "  -c CMD    Execute a shortcut command and exit\n");
    fprintf(stderr, "  -r        Raw JSON mode (pass stdin as-is)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Shortcut commands (-c):\n");
    fprintf(stderr, "  info            Get graph info\n");
    fprintf(stderr, "  cov <kmer>      Query k-mer coverage\n");
    fprintf(stderr, "  find <index>    Find path from index\n");
    fprintf(stderr, "  findkmer <kmer> Find path from k-mer\n");
    fprintf(stderr, "  shutdown        Shutdown server\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Interactive commands (same as shortcuts, or raw JSON)\n");
}

int main(int argc, char *argv[])
{
    const char *host = DEFAULT_HOST;
    int port = DEFAULT_PORT;
    const char *shortcut_cmd = NULL;
    int raw_mode = 0;
    int opt;

    while ((opt = getopt(argc, argv, "h:p:c:r")) != -1) {
        switch (opt) {
        case 'h': host = optarg; break;
        case 'p': port = atoi(optarg); break;
        case 'c': shortcut_cmd = optarg; break;
        case 'r': raw_mode = 1; break;
        default:
            print_usage();
            return 1;
        }
    }

    int fd = connect_to_server(host, port);
    if (fd < 0) {
        fprintf(stderr, "Failed to connect to %s:%d\n", host, port);
        return 1;
    }

    fprintf(stderr, "Connected to DBGPS2-RPC server at %s:%d\n", host, port);

    int rc = 0;

    /* Shortcut command mode */
    if (shortcut_cmd) {
        rc = send_shortcut(fd, shortcut_cmd, 1);
        close(fd);
        return rc < 0 ? 1 : 0;
    }

    /* Raw mode: read JSON from stdin and send */
    if (raw_mode) {
        char line[BUF_SIZE];
        if (fgets(line, sizeof(line), stdin)) {
            /* Strip trailing newline */
            int len = (int)strlen(line);
            while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
                line[--len] = '\0';
            rc = send_request(fd, line, 0);
        }
        close(fd);
        return rc < 0 ? 1 : 0;
    }

    /* Interactive REPL mode */
    printf("DBGPS2-RPC Client - type 'help' for commands, 'quit' to exit\n\n");

    char input[BUF_SIZE];
    int req_id = 1;

    while (1) {
        printf("dbgps2-rpc> ");
        fflush(stdout);

        if (!fgets(input, sizeof(input), stdin)) break;

        /* Strip trailing newline */
        int len = (int)strlen(input);
        while (len > 0 && (input[len-1] == '\n' || input[len-1] == '\r'))
            input[--len] = '\0';

        if (len == 0) continue;
        if (strcmp(input, "quit") == 0 || strcmp(input, "exit") == 0) break;

        if (strcmp(input, "help") == 0) {
            printf("Commands:\n");
            printf("  info            Get graph info\n");
            printf("  cov <kmer>      Query k-mer coverage\n");
            printf("  find <index>    Find path from index\n");
            printf("  findkmer <kmer> Find path from k-mer\n");
            printf("  shutdown        Shutdown server\n");
            printf("  quit / exit     Exit client\n");
            printf("  help            Show this help\n");
            printf("\nOr type a raw JSON-RPC request:\n");
            printf("  {\"jsonrpc\":\"2.0\",\"method\":\"graph/info\",\"id\":1}\n\n");
            continue;
        }

        /* Check if it starts with '{' — raw JSON mode */
        if (input[0] == '{') {
            rc = send_request(fd, input, 1);
        } else {
            /* Shortcut command */
            rc = send_shortcut(fd, input, req_id);
        }

        if (rc < 0) break;
        req_id++;
    }

    close(fd);
    return 0;
}
