# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: default merced check clean realclean

default: build

build: 
	cd ../Src; $(MAKE)
	cp ../Src/calcUpscatterKernel .

check:

clean:
	rm -fr calcUpscatterKernel 

realclean: clean
