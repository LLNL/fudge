# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default check clean realclean doTarget

subDirectories = exponentialIntegral gammaFunctions	angularMomentumCoupling		# erf

default: 

check:
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) doTarget MODE=checke

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
