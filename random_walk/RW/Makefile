
CFLAGS = -std=c99 -g

SRC = ../src/ran.c ../src/config.c ../src/random_walk.c ../src/data.c
SRC_SUM = ../src/ran.c ../src/config.c ../src/summary.c ../src/data.c

random_walk:
	gcc $(SRC) -o rw $(CFLAGS) $(INCLUDES) $(LIBS)

summary:
	gcc $(SRC) -o summary $(CFLAGS) $(INCLUDES) $(LIBS)