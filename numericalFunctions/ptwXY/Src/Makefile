# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(ptwCFLAGS) $(ptwXYCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)
objects = ptwXY_core.o ptwXY_methods.o ptwXY_misc.o ptwXY_functions.o ptwXY_convenient.o ptwXY_integration.o \
    ptwXY_binaryOperators.o ptwXY_unitaryOperators.o ptwXY_interpolation.o

.PHONY: default clean realclean

default: libptwXY.a

libptwXY.a: $(objects)
	$(AR) -cr libptwXY.a *.o

ptwXY_core.o: ptwXY_core.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_core.o ptwXY_core.c

ptwXY_methods.o: ptwXY_methods.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_methods.o ptwXY_methods.c

ptwXY_functions.o: ptwXY_functions.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_functions.o ptwXY_functions.c

ptwXY_misc.o: ptwXY_misc.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_misc.o ptwXY_misc.c

ptwXY_convenient.o: ptwXY_convenient.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_convenient.o ptwXY_convenient.c

ptwXY_integration.o: ptwXY_integration.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_integration.o ptwXY_integration.c

ptwXY_binaryOperators.o: ptwXY_binaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_binaryOperators.o ptwXY_binaryOperators.c

ptwXY_unitaryOperators.o: ptwXY_unitaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_unitaryOperators.o ptwXY_unitaryOperators.c

ptwXY_interpolation.o: ptwXY_interpolation.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_interpolation.o ptwXY_interpolation.c

clean:
	rm -f *.o *.a

realclean: clean
