# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python

.PHONY: default check checke clean realclean

TARGETS = setting.py pop.py domain.py range.py

default:

check:
	for f in $(TARGETS); do $(PYTHON) $$f; done

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:

realclean: clean
