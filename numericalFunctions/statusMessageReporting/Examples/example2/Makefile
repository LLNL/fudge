# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

STATUSMESSAGEREPORTING_CXXFLAGS = -g -O0
local_CXXFLAGS = $(CXXFLAGS) -I./ -I../../Src $(STATUSMESSAGEREPORTING_CXXFLAGS)

.PHONY: default check example clean realclean

default : example

check: default
	echo "******* running example2 *******"
	./example
	echo

example: example.c
	$(CXX) $(local_CXXFLAGS) -o example example.c ../../Src/statusMessageReporting.o

clean:
	rm -f example *.o

realclean: clean
	rm -fr *.dSYM
