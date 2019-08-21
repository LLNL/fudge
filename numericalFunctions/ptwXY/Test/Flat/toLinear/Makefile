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

TARGETS = flatToLinear flatToLinear2

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

flatToLinear: flatToLinear.c
	$(CC) -o flatToLinear flatToLinear.c -lptwXY -lm $(localCFLAGS)

flatToLinear2: flatToLinear2.c
	$(CC) -o flatToLinear2 flatToLinear2.c -lptwXY -lm $(localCFLAGS)

check: default
	for f in $(TARGETS); do $$f -v $(CHECKOPTIONS) > v; python check.py $$f; done

clean:
	rm -f v

realclean: clean
	rm -f $(TARGETS)
