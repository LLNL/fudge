# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTWXY_PATH = ../..
PTWXY_INCLUDE_PATH = $(PTWXY_PATH)/include
PTWXY_LIB_PATH = $(PTWXY_PATH)/lib
PTWXYUTILITIES_LIB_PATH = $(PTWXY_PATH)/Test/Utilities

CDEBUGFLAGS = -g -O0 -Wall
ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTWXY_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -L$(PTWXY_LIB_PATH) \
    -L$(PTWXYUTILITIES_LIB_PATH) -lptwXY -lptwXY_utilities -lnumericalfunctions -lm
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = testBinaryUnion testOverwrite testAddRandomPoints testAddRandomPoints2 testInterpolation_1 testMultiPly2_1 testDivide2_1 testExp_1

.PHONY: default check clean realclean

default:
	$(MAKE) $(TARGETS)

check: default
	for f in $(TARGETS); do python check.py $$f; done

testBinaryUnion: testBinaryUnion.c
	$(CC) -o testBinaryUnion testBinaryUnion.c $(localCFLAGS)

testOverwrite: testOverwrite.c
	$(CC) -o testOverwrite testOverwrite.c $(localCFLAGS)

testAddRandomPoints2: testAddRandomPoints2.c
	$(CC) -o testAddRandomPoints2 testAddRandomPoints2.c $(localCFLAGS)

testAddRandomPoints: testAddRandomPoints.c
	$(CC) -o testAddRandomPoints testAddRandomPoints.c $(localCFLAGS)

testInterpolation_1: testInterpolation_1.c
	$(CC) -o testInterpolation_1 testInterpolation_1.c $(localCFLAGS)

testMultiPly2_1: testMultiPly2_1.c
	$(CC) -o testMultiPly2_1 testMultiPly2_1.c $(localCFLAGS)

testDivide2_1: testDivide2_1.c
	$(CC) -L. -o testDivide2_1 testDivide2_1.c $(localCFLAGS)

testExp_1: testExp_1.c
	$(CC) -L. -o testExp_1 testExp_1.c $(localCFLAGS)

clean:
	rm -f *.dat t_*

realclean: clean
	rm -f $(TARGETS)
