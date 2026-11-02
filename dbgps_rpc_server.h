#ifndef DBGPS_RPC_SERVER_H
#define DBGPS_RPC_SERVER_H

#include "dbgps_core.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    dbgps_graph_t *graph;
    int port;
    int max_clients;
    volatile int running;
} dbgps_rpc_server_t;

/**
 * Start the DBGPS2-RPC server. This function blocks until the server
 * is stopped (via a "server/shutdown" JSON-RPC command or SIGINT).
 *
 * Returns 0 on clean shutdown, -1 on error.
 */
int dbgps_rpc_server_start(dbgps_rpc_server_t *srv);

/**
 * Signal the server to stop. Can be called from a signal handler
 * or another thread.
 */
void dbgps_rpc_server_stop(dbgps_rpc_server_t *srv);

#ifdef __cplusplus
}
#endif

#endif /* DBGPS_RPC_SERVER_H */
