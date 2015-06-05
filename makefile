CC=gcc
CFLAGS= -g

calc: main.c bessel_mod.c integrate.c
	$(CC) -o calc bessel_mod.c  main.c integrate.c -lm -lgsl -lgslcblas 
