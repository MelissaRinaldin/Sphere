CC=gcc

all:
	$(CC) -o sphere sphere.c -O3 -lm
