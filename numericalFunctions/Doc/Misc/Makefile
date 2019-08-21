# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../../
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

CC = gcc
CDEBUGFLAGS = -g -O0 -Wall
ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = infill_f3_f4

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

infill_f3_f4: infill_f3_f4.c
	$(CC) -o infill_f3_f4 infill_f3_f4.c -lm $(localCFLAGS)

check:

clean:
	rm -f $(TARGETS)

realclean: clean
