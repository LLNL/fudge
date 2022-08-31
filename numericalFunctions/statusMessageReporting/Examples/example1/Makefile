# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

STATUSMESSAGEREPORTING_CXXFLAGS = -g -O0
local_CXXFLAGS = $(CXXFLAGS) -I./ -I../../Src $(STATUSMESSAGEREPORTING_CXXFLAGS)

.PHONY: all check clean realclean

default : example

check: default
	echo "******* running example1 *******"
	./example
	echo
	echo "******* re-running example1 with argument *******"
	./example 1
	echo

example: example.c checkMathSMR.o mathSMR.o
	$(CXX) $(local_CXXFLAGS) -o example example.c checkMathSMR.c mathSMR.o ../../Src/statusMessageReporting.o -lm

checkMathSMR.o: checkMathSMR.c
	$(CXX) $(local_CXXFLAGS) -c -o checkMathSMR.o checkMathSMR.c

mathSMR.o: mathSMR.c
	$(CXX) $(local_CXXFLAGS) -c -o mathSMR.o mathSMR.c

clean:
	rm -f example *.o

realclean: clean
	rm -fr *.dSYM
