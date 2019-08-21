# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh
TARGETS = gammaTest incompleteGammaTest

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib


CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I./

.PHONY: default clean realclean m gammaMathematica.m incompleteGammaMathematica.m

default: $(TARGETS)

check: default
	./gammaTest $(CHECKOPTIONS)
	./incompleteGammaTest $(CHECKOPTIONS)

m: gammaMathematica.m incompleteGammaMathematica.m

gammaTest: gammaTest.c
	$(CC) $(localCFLAGS) -o gammaTest gammaTest.c -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions -lm

incompleteGammaTest: incompleteGammaTest.c
	$(CC) $(localCFLAGS) -o incompleteGammaTest incompleteGammaTest.c -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions -lm

gammaMathematica.m:
	$(CC) $(CFLAGS) -o gammaMathematica gammaMathematica.c
	./gammaMathematica > gammaMathematica.m

incompleteGammaMathematica.m:
	$(CC) $(CFLAGS) -o incompleteGammaMathematica incompleteGammaMathematica.c
	./incompleteGammaMathematica > incompleteGammaMathematica.m

clean:
	rm -f gammaTest gammaMathematica incompleteGammaTest incompleteGammaMathematica *.m

realclean: clean
