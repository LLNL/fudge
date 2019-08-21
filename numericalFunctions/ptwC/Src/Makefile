# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../..

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS)
objects = nf_utilities.o

.PHONY: default clean realclean

default: libnf_utilities.a

libnf_utilities.a: $(objects)
	$(AR) -cr libnf_utilities.a *.o

nf_utilities.o: nf_utilities.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_utilities.o nf_utilities.c

clean:
	rm -f *.o *.a

realclean: clean
