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

TARGETS = xArray clip trim dullEdges mutualify mod areDomainsMutual mergeClosePoints copy

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

xArray: xArray.c
	$(CC) -o xArray xArray.c -lptwXY -lm $(localCFLAGS)

clip: clip.c
	$(CC) -o clip clip.c -lptwXY -lm $(localCFLAGS)

trim: trim.c
	$(CC) -o trim trim.c -lptwXY -lm $(localCFLAGS)

dullEdges: dullEdges.c
	$(CC) -o dullEdges dullEdges.c -lptwXY -lm $(localCFLAGS)

mutualify: mutualify.c
	$(CC) -o mutualify mutualify.c -lptwXY -lm $(localCFLAGS)

mod: mod.c
	$(CC) -o mod mod.c -lptwXY -lm $(localCFLAGS)

areDomainsMutual: areDomainsMutual.c
	$(CC) -o areDomainsMutual areDomainsMutual.c -lptwXY -lm $(localCFLAGS)

mergeClosePoints: mergeClosePoints.c
	$(CC) -o mergeClosePoints mergeClosePoints.c -lptwXY -lm $(localCFLAGS)

copy: copy.c
	$(CC) -o copy copy.c -lptwXY -lm $(localCFLAGS)

check: default
	for f in $(TARGETS); do $$f $(CHECKOPTIONS); done

clean:
	rm -f v

realclean: clean
	rm -f $(TARGETS)
