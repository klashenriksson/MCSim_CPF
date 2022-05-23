
CFLAGS = -std=c99 -O3

SRC = ../src/ran.c ../src/config.c ../src/percolation.c ../src/data.c ../src/int_queue.c
SRC_SUM = ../src/ran.c ../src/config.c ../src/summary.c ../src/data.c

random_walk:
	gcc $(SRC) -o perc $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc $(SRC_SUM) -o summary $(CFLAGS) $(INCLUDES) $(LIBS)