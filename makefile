CC=gcc
CFLAGS=-Wall -lm

calc: main.c bessel_mod.c 
	$(CC) -o calc bessel_mod.c  main.c -lm 
