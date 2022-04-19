CC = gcc
CPPFLAGS = -DCLU
CFLAGS  = -g -O4 -std=c99
LOADLIBES = -lm
VPATH = ../src

INCLUDES = ran.h

OBJS = ising.o ran.o config.o corr.o int_queue.o

ising: ${OBJS}

summary: summary.o

${OBJS}: ${INCLUDES} Makefile
