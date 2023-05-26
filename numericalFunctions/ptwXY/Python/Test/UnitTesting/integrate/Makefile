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
	$(PYTHON) groupTestAll1.py
	$(PYTHON) integrationXY.py
	$(PYTHON) integrationWithWeight_x.py
	$(PYTHON) integrationWithWeight_sqrt_x.py
	$(PYTHON) integration2.py
	$(PYTHON) runningIntegration.py

clean:

realclean: clean
