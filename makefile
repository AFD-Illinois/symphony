CC=gcc -std=c99
CFLAGS= -Wall -Ofast

symphony: main.c bessel_mod.c dec.h params.h calc.c fits.c
	$(CC) -o symphony params.h dec.h bessel_mod.c calc.c main.c fits.c -lm -lgsl -lgslcblas 
