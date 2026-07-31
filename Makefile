CC=gcc --std=c99
CFLAGS=-O2

.PHONY: all

all: finite_vol

finite_vol: finite_vol.c
	$(CC) $(CFLAGS) $< -o $@ -lm 


