# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

subDirectories = Utilities initial UnitTesting Flat

.PHONY: default utilities check checke clean realclean doTarget

default:
	$(MAKE) doTarget MODE=''

utilities:
	cd Utilities; $(MAKE) utilities

check: default
	-$(MAKE) doTarget MODE=check

checke:
	$(MAKE) doTarget MODE=check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
