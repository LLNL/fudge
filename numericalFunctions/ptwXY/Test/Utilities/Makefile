# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: check clean realclean

CC = gcc
CFLAGS = -g -O0 -Wall
ARFLAGS = rc

NUMERICALFUNCTIONS_PATH = ../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

PTWXY_PATH = ../..
PTWXY_INCLUDE_PATH = $(PTWXY_PATH)/include

localCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTWXY_INCLUDE_PATH)

libptwXY_utilities.a: ptwXY_utilities.o
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

ptwXY_utilities.o: ptwXY_utilities.c ptwXY_utilities.h
	$(CC) $(CFLAGS) -c -o ptwXY_utilities.o ptwXY_utilities.c $(localCFLAGS)

check: libptwXY_utilities.a

clean:
	rm -f *.a *.o

realclean: clean
