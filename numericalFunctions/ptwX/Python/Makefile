# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default clean realclean doTarget

subDirectories = Src Test

default:
	cd Src; $(MAKE)

check:
	cd Test; $(MAKE) check

checke:
	cd Test; $(MAKE) check CHECKOPTIONS=-e

clean:
	rm -rf ptwXY ptwC ptwX
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf ptwXY ptwC ptwX
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
