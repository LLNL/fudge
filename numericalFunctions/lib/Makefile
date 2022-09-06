# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default libnumericalfunctions.a check clean realclean

default: libnumericalfunctions.a

libnumericalfunctions.a:
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
	rm -fr *.pyc __.SYMDEF* __pycache__
