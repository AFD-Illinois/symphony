CC=gcc
CFLAGS= -Ofast -Wall

symphony: main.c bessel_mod.c dec.h params.h
	$(CC) -o symphony params.h dec.h bessel_mod.c main.c -lm -lgsl -lgslcblas 
