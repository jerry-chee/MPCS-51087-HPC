# Jerry Chee
# MPCS 51807 hw2 
# Makefile

CC := mpicc
CFLAGS := -std=c99 -Wall -g -O0
OBJ := mandelbrot.o static_decomp.o dynamic_decomp.o complex.o
SRC := mandelbrot.c static_decomp.c dynamic_decomp.c complex.c

mandelbrot: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o mandelbrot

mandelbrot.o: $(SRC) decomp.h complex.h
	$(CC) $(CFLAGS) -c mandelbrot.c

static_decomp.o: static_decomp.c decomp.h
	$(CC) $(CFLAGS) -c static_decomp.c

dynamic_decomp.o: dynamic_decomp.c static_decomp.c decomp.h
	$(CC) $(CFLAGS) -c dynamic_decomp.c

complex.o: complex.c complex.h
	$(CC) $(CFLAGS) -c complex.c

latency_bandwidth: latency_bandwidth.c
	$(CC) $(CFLAGS) latency_bandwidth.c -o latency_bandwidth
