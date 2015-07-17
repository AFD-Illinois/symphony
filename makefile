CC=gcc
CFLAGS= -Ofast -Wall

symphony: main.c bessel_mod.c dec.h params.h calc.c
	$(CC) -o symphony params.h dec.h bessel_mod.c calc.c main.c -lm -lgsl -lgslcblas 
