#
# Makefile for raytr
# L. Benkevitch March 10 2011
#
CFLAGS = -g -Wall
#CPPFLAGS =
LDLIBS = -lm

CC = gcc

write_rwimage: write_rwimage.o rwimage.o rwfun.o beam_path.o


install: write_rwimage
 
clean:
	 rm -f *.o 
