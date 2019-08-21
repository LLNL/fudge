# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh
TARGETS = exponentialIntegralTest

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib


CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I./

.PHONY: default clean realclean

default: $(TARGETS)

check: default
	./exponentialIntegralTest $(CHECKOPTIONS)

exponentialIntegralTest: exponentialIntegralTest.c
	$(CC) $(localCFLAGS) -o exponentialIntegralTest exponentialIntegralTest.c -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions -lm

exponentialIntegralMathematica:
	$(CC) $(CFLAGS) -o exponentialIntegralMathematica exponentialIntegralMathematica.c
	./exponentialIntegralMathematica > exponentialIntegralMathematica.m

clean:
	rm -f exponentialIntegralTest exponentialIntegralMathematica exponentialIntegralMathematica.m

realclean: clean
