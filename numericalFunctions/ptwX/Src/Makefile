# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

PTWC_PATH = ../../ptwC
PTWC_INCLUDE_PATH = $(PTWC_PATH)/include

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(ptwCFLAGS) $(ptwXCFLAGS) -I$(PTWC_INCLUDE_PATH)
objects = ptwX_core.o

.PHONY: default check clean realclean

default: libptwX.a

check:

libptwX.a: $(objects)
	$(AR) -cr libptwX.a *.o

ptwX_core.o: ptwX_core.c ptwX.h
	$(CC) $(localCFLAGS) -c -o ptwX_core.o ptwX_core.c

clean:
	rm -f *.o *.a

realclean: clean
