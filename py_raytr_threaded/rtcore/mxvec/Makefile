#
# Makefile for libmxv.a, a small matrix-vector library
# used in the raytrace Python extension module
# Created March 15, 2011 by Leonid Benkevitch. 
#
CFLAGS = -g -fPIC -Wall
CPPFLAGS = -I..
LDLIBS = -lm

CC = gcc

MXVOBJS=\
cross_product.o dot_product.o linsolve.o linsolve2.o maxval.o minval.o \
minv.o minv2.o mmul.o mvmul.o sum.o sum_squares.o print1d.o print2d.o \
vmagn.o v2magn.o v3magn.o logspace.o log10space.o calc_rpos.o calc_rdir.o \
mcalc_rdir.o mcalc_rsph.o

libmxv.a: $(MXVOBJS)
	ar cvr libmxv.a $(MXVOBJS)
	ranlib libmxv.a

install: libmxv.a
	cp libmxv.a   ../lib

clean:
	 rm -f *.o libmxv.a

purge:
	@echo 'Cleaning remporary emacs *~ and vi #*# files'; \
	rm -f `find ./ -name '*~'` && rm -f `find ./ -name '#*#'`; 
	rm -f *.o *.a