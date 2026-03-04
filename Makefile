CC=mpicc
CC_EXAMPLE=gcc
CFLAGS=-O2 -Wall -Wextra -std=c11 -Iinclude
LDFLAGS=-lm

SRC=$(wildcard src/*.c)
OBJ=$(SRC:.c=.o)

all: stokes

.PHONY: all clean fluid3d

stokes: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

fluid3d: examples/fluid3d_semi_lagrangian.c
	$(CC_EXAMPLE) -O2 -Wall -Wextra -std=c11 $< -o $@ -lm

clean:
	rm -f src/*.o stokes fluid3d