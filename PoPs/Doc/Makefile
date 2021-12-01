# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/bash

.PHONY: default doc clean realclean

default: doc

doc:
	cd sphinx; $(MAKE) html
	cp -r sphinx/build/html .

check:

clean:
	cd sphinx; $(MAKE) clean
	rm -rf html

realclean: clean
