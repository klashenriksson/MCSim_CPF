
INCLUDES = -I../g2/include -L../g2
LIBS = -lm
CFLAGS = -std=c99 -DCORR -O3

SRC = ../src/ran.c ../src/config.c ../src/ising.c ../src/corr.c

ising:
	gcc $(SRC) -o ising $(CFLAGS) $(INCLUDES) $(LIBS)