# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python3

localCFLAGS = $(CFLAGS) -I../include

.PHONY: all check clean realclean

default :

check: check1 check2 check3
	./check1
	./check2 > check2.out
	diff check2.out Data/check2.out > /dev/null; if [ $$? != 0 ]; then echo "FAILED: check2"; fi
	./check3 > check3.out
	$(PYTHON) check3.py
	diff check3.out Data/check3.out > /dev/null; if [ $$? != 0 ]; then echo "FAILED: check3"; fi

check1: check1.c
	$(CC) $(localCFLAGS) -o check1 check1.c -L../lib -lstatusMessageReporting

check2: check2.c
	$(CC) $(localCFLAGS) -o check2 check2.c -L../lib -lstatusMessageReporting

check3: check3.c
	$(CC) $(localCFLAGS) -o check3 check3.c -L../lib -lstatusMessageReporting

clean:
	rm -f check1 check2 check2.out check3 check3.out

realclean: clean
	rm -fr *.dSYM
