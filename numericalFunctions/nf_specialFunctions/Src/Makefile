# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)

objects = nf_exponentialIntegral.o nf_gammaFunctions.o nf_incompleteGammaFunctions.o nf_polevl.o

.PHONY: default clean realclean

default: $(objects)

nf_exponentialIntegral.o: nf_exponentialIntegral.c
	$(CC) $(localCFLAGS) -c -o nf_exponentialIntegral.o nf_exponentialIntegral.c

nf_gammaFunctions.o: nf_gammaFunctions.c
	$(CC) $(localCFLAGS) -c -o nf_gammaFunctions.o nf_gammaFunctions.c

nf_incompleteGammaFunctions.o: nf_incompleteGammaFunctions.c
	$(CC) $(localCFLAGS) -c -o nf_incompleteGammaFunctions.o nf_incompleteGammaFunctions.c

nf_polevl.o: nf_polevl.c
	$(CC) $(localCFLAGS) -c -o nf_polevl.o nf_polevl.c

clean:
	rm -f *.o *.a

realclean: clean
