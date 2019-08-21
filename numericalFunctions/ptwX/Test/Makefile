# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check clean realclean doTarget

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTW1D_PATH = ..
PTW1D_INCLUDE_PATH = $(PTW1D_PATH)/include
PTW1D_LIB_PATH = $(PTW1D_PATH)/lib
PTW1D_ARCHIVE_FILE = $(PTW1D_LIB_PATH)/libptwX.a

CDEBUGFLAGS = -g -O0 -Wall
ptwXCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTW1D_INCLUDE_PATH) -L$(PTW1D_LIB_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lptwX -lnumericalfunctions -lm
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXCFLAGS)

subDirectories = 

default: testCreate

testCreate: testCreate.c $(PTWXY_ARCHIVE_FILE)
	$(CC) -o testCreate testCreate.c $(localCFLAGS)

check: default
	./testCreate $(CHECKOPTIONS)

clean:
	rm -f testCreate
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
