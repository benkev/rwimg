#
# Makefile for raytr
# L. Benkevitch March 10 2011
#
CFLAGS = -g -Wall
#CPPFLAGS =
LDFLAGS = -L.
LDLIBS = -lmxv -lm

CC = gcc

#c.o.:
#	 gcc $(CFLAGS) -c $<

SOURCES=\
write_rwimage.c rwimage.c beam_path.c plasma_density_inv_squares.c

OBJECTS=\
write_rwimage.o rwimage.o beam_path.o plasma_density_inv_squares.o

MXVOBJS=\
cross_product.o dot_product.o linsolve.o linsolve2.o maxval.o minval.o \
minv.o minv2.o mmul.o mvmul.o sum.o sum_squares.o print1d.o  print2d.o

raytr: libmxv.a $(OBJECTS) 
	 gcc $(OBJECTS) -L. $(LDLIBS) -o raytr

install: raytr

libmxv.a: $(MXVOBJS)
	ar cvr libmxv.a $(MXVOBJS)
	ranlib libmxv.a

clean:
	 rm -f *.o 
