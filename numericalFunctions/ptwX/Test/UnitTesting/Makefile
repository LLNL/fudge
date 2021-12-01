# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default check clean realclean doTarget

subDirectories = creating binary

default:
	$(MAKE) MODE=default doTarget

check: default
	$(MAKE) $(CHECKOPTIONS) MODE=check doTarget

clean:
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
