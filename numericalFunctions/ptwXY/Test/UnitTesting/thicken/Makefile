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

TARGETS = thicken

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

check: default
	./thicken $(CHECKOPTIONS)

thicken: thicken.c
	$(CC) -o thicken thicken.c -lptwXY -lm $(localCFLAGS)

clean:
	rm -f $(TARGETS) v

realclean: clean
