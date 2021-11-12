# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default utilities check checke checkPython clean realclean doTarget

subDirectories = Python Src Test lib include

default:
	cd lib; $(MAKE)
	cd include; $(MAKE)

utilities:
	cd Test; $(MAKE) utilities

check: default
	cd Test; $(MAKE) check

checke:
	$(MAKE) check CHECKOPTIONS=-e

checkPython:
	cd Python; $(MAKE) check

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
