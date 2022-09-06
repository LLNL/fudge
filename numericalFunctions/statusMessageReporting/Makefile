# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default lib include check clean realclean doSubDirs

subDirs = Src Test lib include Doc _Src Examples

default: include lib

lib:
	cd lib; $(MAKE)

include:
	cd include; $(MAKE)

check: default
	cd Test; $(MAKE) check

collectSources:
	cd _Src; $(MAKE) collectSources

clean:
	$(MAKE) doTarget=clean doSubDirs

realclean:
	$(MAKE) doTarget=realclean doSubDirs

doSubDirs:
	pwDir=`pwd`; \
	for subDir in $(subDirs); do cd $$pwDir/$$subDir; $(MAKE) $(doTarget); done
