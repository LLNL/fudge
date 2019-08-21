# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default libnumericalfunctions.a check clean realclean

default: libnumericalfunctions.a

libnumericalfunctions.a: $(objects)
	rm -f *.o
	$(AR) -x $(NFUPDATING_PATH)/lib/*.a
	ARFLAGS=-cr; \
	if [ -e libnumericalfunctions.a ]; then ARFLAGS=-r; fi; \
	$(AR) $$ARFLAGS libnumericalfunctions.a *.o;
	rm -f *.o

check:

clean:
	rm -f *.o *.a *.so

realclean: clean
	rm -f *.pyc
