# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default check checke clean realclean doTarget

subDirectories = thicken integrate Others convolution thinning binaryMath interpolation toFunction inverse

default: 
	$(MAKE) doTarget

check:
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) doTarget MODE=check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
