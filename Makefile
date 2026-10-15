CFLAGS=-g -Wall -O3
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=DBGPS2 DBGPS3 dbgps3-client

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean version test

all:$(PROG)

# Keep test executables outside the checkout. Use `make test asan=1` for ASan.
test:
	@set -eu; test_dir=$$(mktemp -d); trap 'rm -rf "$$test_dir"' EXIT; \
	for source in DBGPS2.c dbgps3_core.c; do \
		defs=; if [ "$$source" = dbgps3_core.c ]; then defs=-DTEST_DBGPS3; fi; \
		$(CC) $(CFLAGS) -I. $$defs -DTEST_SOURCE=\"$$source\" -o "$$test_dir/cov" tests/test_kmer_cov.c dbgps_codec.c kthread.c $(LIBS) -lpthread; \
		"$$test_dir/cov"; \
	done; \
	for defs in '' '-DTEST_DBGPS3'; do \
		$(CC) $(CFLAGS) -I. $$defs -o "$$test_dir/paths" tests/test_paths.c dbgps_codec.c kthread.c $(LIBS) -lpthread; \
		"$$test_dir/paths"; \
	done; \
	$(CC) $(CFLAGS) -I. -o "$$test_dir/codec" tests/test_codec.c dbgps_codec.c; \
	"$$test_dir/codec" "$$test_dir/generated.map"; \
	$(CC) $(CFLAGS) -o "$$test_dir/DBGPS2" DBGPS2.c dbgps_codec.c kthread.c $(LIBS) -lpthread; \
	$(CC) $(CFLAGS) -o "$$test_dir/DBGPS3" DBGPS3.c dbgps3_core.c dbgps3_server.c minjson.c kthread.c $(LIBS) -lpthread; \
	python3 tests/test_cli.py "$$test_dir/DBGPS2" "$$test_dir/DBGPS3"; \
	python3 tests/test_mapping_cli.py "$$test_dir/DBGPS2"

# DBGPS2: Read-Length constrained version with dynamic coverage ratio support
# Previously DBGPS-dy2-rl
# Features:
# - Skip option for initial ratio checks
# - Output prefix option for log and strands files
# - Ratio range iteration support
# - Improved path finding algorithm
DBGPS2:DBGPS2.c dbgps_codec.c dbgps_codec.h khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS2.c dbgps_codec.c kthread.c $(LIBS) -lpthread

# DBGPS3: Modular version with CLI + server mode
# Features:
# - Refactored into dbgps3_core (graph engine) + dbgps3_server (JSON-RPC server)
# - Server mode via -S flag with JSON-RPC 2.0 over TCP
# - Thread-per-connection concurrency model
# - minjson: lightweight JSON parser/builder for JSON-RPC protocol
DBGPS3:DBGPS3.c dbgps3_core.c dbgps3_server.c minjson.c khashl.h ketopt.h kseq.h kthread.h dbgps3_core.h dbgps3_server.h minjson.h
	$(CC) $(CFLAGS) -o $@ DBGPS3.c dbgps3_core.c dbgps3_server.c minjson.c kthread.c $(LIBS) -lpthread

# dbgps3-client: Interactive CLI client for DBGPS3 server
dbgps3-client:dbgps3_client.c minjson.c minjson.h
	$(CC) $(CFLAGS) -o $@ dbgps3_client.c minjson.c $(LIBS)

clean:
	rm -fr *.dSYM $(PROG)

# Build a specific version: make VERSION=DBGPS2
VERSION?=all
version:
	$(MAKE) $(VERSION)

