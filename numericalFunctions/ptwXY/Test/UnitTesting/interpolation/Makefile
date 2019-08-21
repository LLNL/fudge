# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTWXY_PATH = ../../..
PTWXY_INCLUDE_PATH = $(PTWXY_PATH)/include
PTWXY_LIB_PATH = $(PTWXY_PATH)/lib
PTWXYUTILITIES_INCLUDE_PATH = $(PTWXY_PATH)/Test/Utilities
PTWXYUTILITIES_LIB_PATH = $(PTWXY_PATH)/Test/Utilities

CC = gcc
CDEBUGFLAGS = -g -O0 -Wall
ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTWXY_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -L$(PTWXY_LIB_PATH) \
    -I$(PTWXYUTILITIES_INCLUDE_PATH) -L$(PTWXYUTILITIES_LIB_PATH) -lptwXY_utilities -lnumericalfunctions
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = toLinearLinear unitbase unitbase2

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

toLinearLinear: toLinearLinear.c
	$(CC) -o toLinearLinear toLinearLinear.c -lptwXY -lm $(localCFLAGS)

unitbase: unitbase.c
	$(CC) -o unitbase unitbase.c -lptwXY -lm $(localCFLAGS)

unitbase2: unitbase2.c
	$(CC) -o unitbase2 unitbase2.c -lptwXY -lm $(localCFLAGS)

check: default
	./toLinearLinear $(CHECKOPTIONS)
	./toLinearLinear $(CHECKOPTIONS) -xlog
	./toLinearLinear $(CHECKOPTIONS) -ylog
	./toLinearLinear $(CHECKOPTIONS) -xlog -ylog
	./unitbase $(CHECKOPTIONS)
	./unitbase2 $(CHECKOPTIONS)

clean:
	rm -f *.dat v

realclean: clean
	rm -f $(TARGETS)
