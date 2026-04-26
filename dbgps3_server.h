#ifndef DBGPS3_SERVER_H
#define DBGPS3_SERVER_H

#include "dbgps3_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    dbgps3_graph_t *graph;
    int port;
    int max_clients;
    volatile int running;
} dbgps3_server_t;

/**
 * Start the DBGPS3 server. This function blocks until the server
 * is stopped (via a "server/shutdown" JSON-RPC command or SIGINT).
 *
 * Returns 0 on clean shutdown, -1 on error.
 */
int dbgps3_server_start(dbgps3_server_t *srv);

/**
 * Signal the server to stop. Can be called from a signal handler
 * or another thread.
 */
void dbgps3_server_stop(dbgps3_server_t *srv);

#ifdef __cplusplus
}
#endif

#endif /* DBGPS3_SERVER_H */
