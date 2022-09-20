# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python

.PHONY: default check clean realclean

default:

check:
	$(PYTHON) thin.py
	$(PYTHON) thinDomain.py

clean:

realclean: clean
