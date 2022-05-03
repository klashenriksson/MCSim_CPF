
INCLUDES = -I../g2/include -L../g2
LIBS = -lm
CFLAGS = -std=c99 -DCORR -DCLU -O3

SRC = ../src/ran.c ../src/config.c ../src/ising.c ../src/corr.c ../src/int_queue.c ../src/data.c

SRC_SUMMARY = ../src/data.c ../src/summary.c

ising:
	gcc $(SRC) -o ising $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc $(SRC_SUMMARY) -o summary $(CFLAGS) $(INCLUDES) $(LIBS)