# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh
PYTHON = python3

.PHONY: default doc check doctest clean realclean

default:

doc:
	cd Doc; $(MAKE)

check:
	cd Check; $(MAKE) check PYTHON=$(PYTHON)

doctest:
	$(PYTHON) -m doctest -v PQU.py

clean:
	cd Doc; $(MAKE) clean
	cd Check; $(MAKE) clean

realclean:
	rm -rf *.pyc
	rm -rf __pycache__
	cd Doc; $(MAKE) realclean
	cd Check; $(MAKE) realclean
