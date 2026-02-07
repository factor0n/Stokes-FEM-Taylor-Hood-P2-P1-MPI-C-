CC=mpicc
CFLAGS=-O2 -Wall -Wextra -std=c11 -Iinclude
LDFLAGS=-lm

SRC=$(wildcard src/*.c)
OBJ=$(SRC:.c=.o)

all: stokes

stokes: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

clean:
	rm -f src/*.o stokes
