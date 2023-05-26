# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python3

.PHONY: default check checke clean realclean

TARGETS = clip.py trim.py xArray.py dullEdges.py mutualify.py mod.py areDomainsMutual.py reallocate.py mergeClosePoints.py \
		copy.py lowerIndexBoundingX.py

default:

check:
	for f in $(TARGETS); do $(PYTHON) $$f; done

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:

realclean: clean
