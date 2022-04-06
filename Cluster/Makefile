CC = gcc
CPPFLAGS = -DCLU
CFLAGS  = -g -O4
LOADLIBES = -lm
VPATH = ../src

INCLUDES = ran.h

OBJS = ising.o ran.o config.o

ising: ${OBJS}

summary: summary.o

${OBJS}: ${INCLUDES} Makefile
