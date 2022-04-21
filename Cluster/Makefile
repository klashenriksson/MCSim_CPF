
INCLUDES = -I../g2/include -L../g2
LIBS = -lm -lg2 -lgd
CFLAGS = -std=c99 -DCORR -DDRAW_GRAPHICS -DCLU -O3

SRC = ../src/ran.c ../src/config.c ../src/ising.c ../src/corr.c ../src/int_queue.c

ising:
	gcc $(SRC) -o ising $(CFLAGS) $(INCLUDES) $(LIBS)