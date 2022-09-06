# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

subDirs = Src bin

.PHONY: default bin clean realclean check doc doSubDirs

default: bin

bin:
	cd bin; $(MAKE)

clean:
	$(MAKE) doTarget=clean doSubDirs

realclean:
	$(MAKE) doTarget=realclean doSubDirs

check:

doSubDirs:
	SAVED_PWD=`pwd`; \
	for subDir in $(subDirs); do cd $$SAVED_PWD/$$subDir; $(MAKE) $(doTarget); done
