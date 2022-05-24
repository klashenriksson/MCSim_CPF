
CFLAGS = -std=c99 -O3 -DSURV_BIAS

SRC = ../src/ran.c ../src/config.c ../src/random_walk.c ../src/data.c
SRC_SUM = ../src/ran.c ../src/config.c ../src/summary.c ../src/data.c

random_walk:
	gcc $(SRC) -o rw $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc $(SRC_SUM) -o summary $(CFLAGS) $(INCLUDES) $(LIBS)