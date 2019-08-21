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

TARGETS = groupTestAll1 integration integrationXY integrationWithWeight_x integrationWithWeight_sqrt_x

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

groupTestAll1: groupTestAll1.c
	$(CC) -o groupTestAll1 groupTestAll1.c -lptwXY -lm $(localCFLAGS)

integration: integration.c
	$(CC) -o integration integration.c -lptwXY -lm $(localCFLAGS)

integrationXY: integrationXY.c
	$(CC) -o integrationXY integrationXY.c -lptwXY -lm $(localCFLAGS)

integrationWithWeight_x: integrationWithWeight_x.c
	$(CC) -o integrationWithWeight_x integrationWithWeight_x.c -lptwXY -lm $(localCFLAGS)

integrationWithWeight_sqrt_x: integrationWithWeight_sqrt_x.c
	$(CC) -o integrationWithWeight_sqrt_x integrationWithWeight_sqrt_x.c -lptwXY -lm $(localCFLAGS)

check: default
	for f in $(TARGETS); do $$f $(CHECKOPTIONS); done

clean:
	rm -f *.dat v

realclean: clean
	rm -f $(TARGETS)
