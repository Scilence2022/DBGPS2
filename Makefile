CFLAGS=-g -Wall -O3
CXXFLAGS=$(CFLAGS) -std=c++11
LIBS=-lz
PROG=DBGPS DBGPS-dy DBGPS-dy2 DBGPS-greedy-dy2 DBGPS-greedy-paths DBGPS-ft DBGPS-SmKnd DBGPS-analyzer DBGPS-pc DBGPS-links DBGPS-ps DBGPS-mt  DBGPS-greedy DBGPS-dy2-rl

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.PHONY:all clean

all:$(PROG)

DBGPS:DBGPS.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS.c kthread.c $(LIBS) -lpthread

DBGPS-dy:DBGPS-dy.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-dy.c kthread.c $(LIBS) -lpthread

DBGPS-dy2:DBGPS-dy2.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-dy2.c kthread.c $(LIBS) -lpthread

DBGPS-greedy:DBGPS-greedy.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-greedy.c kthread.c $(LIBS) -lpthread

DBGPS-greedy-dy2:DBGPS-greedy-dy2.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-greedy-dy2.c kthread.c $(LIBS) -lpthread

DBGPS-greedy-paths:DBGPS-greedy-paths.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-greedy-paths.c kthread.c $(LIBS) -lpthread

DBGPS-ft:DBGPS-ft.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-ft.c kthread.c $(LIBS) -lpthread

DBGPS-SmKnd:DBGPS-SmKnd.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-SmKnd.c kthread.c $(LIBS) -lpthread

DBGPS-analyzer:DBGPS-analyzer.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-analyzer.c kthread.c $(LIBS) -lpthread

DBGPS-pc:DBGPS-pc.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-pc.c kthread.c $(LIBS) -lpthread

DBGPS-links:DBGPS-links.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-links.c kthread.c $(LIBS) -lpthread

DBGPS-ps:DBGPS-ps.c khashl.h ketopt.h kseq.h kthread.h fpthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-ps.c kthread.c fpthread.c $(LIBS) -lpthread

DBGPS-mt:DBGPS-mt.c khashl.h ketopt.h kseq.h kthread.h fpthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-mt.c kthread.c fpthread.c $(LIBS) -lpthread

DBGPS-dy2-rl:DBGPS-dy2-rl.c khashl.h ketopt.h kseq.h kthread.h
	$(CC) $(CFLAGS) -o $@ DBGPS-dy2-rl.c kthread.c $(LIBS) -lpthread

clean:
	rm -fr *.dSYM $(PROG)

