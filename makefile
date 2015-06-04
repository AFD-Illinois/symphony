CC=gcc
CFLAGS=-Wall 

calc: main.c bessel_mod.c 
	$(CC) -o calc bessel_mod.c  main.c -lm -lgsl -lgslcblas
