# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python3

.PHONY: default check clean realclean

default:

check:
	$(PYTHON) getSlice.py
	$(PYTHON) setSlice.py
	$(PYTHON) xSlice.py

checke:
	$(MAKE) CHECKOPTIONS="-v" check

clean:

realclean: clean
