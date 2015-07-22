CC=gcc -std=c99
CFLAGS= -Wall -Ofast

symphony: main.c bessel_mod.c symphony.h params.h calc.c fits.c integrate.c
	$(CC) -o symphony params.h symphony.h bessel_mod.c calc.c main.c fits.c integrate.c -lm -lgsl -lgslcblas 
