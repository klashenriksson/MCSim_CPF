CC = gcc
CFLAGS  = -g -O4
LOADLIBES = -lm
VPATH = ../src

INCLUDES = ran.h

OBJS = ising.o config.o ran.o

ising: ${OBJS}

summary.o: ising.h

summary: summary.o

gtot: gtot.o

${OBJS}: ${INCLUDES} Makefile
