# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(nf_LegendreCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)
objects = nf_Legendre.o nf_Legendre_GaussianQuadrature.o

.PHONY: default clean realclean

default: libnf_Legendre.a

libnf_Legendre.a: $(objects)
	$(AR) -cr libnf_Legendre.a *.o

nf_Legendre.o: nf_Legendre.c nf_Legendre.h
	$(CC) $(localCFLAGS) -c -o nf_Legendre.o nf_Legendre.c

nf_Legendre_GaussianQuadrature.o: nf_Legendre_GaussianQuadrature.c nf_Legendre.h
	$(CC) $(localCFLAGS) -c -o nf_Legendre_GaussianQuadrature.o nf_Legendre_GaussianQuadrature.c

clean:
	rm -f *.o *.a

realclean: clean
