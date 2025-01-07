# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: default merced check clean realclean

default: merced upscatter

merced:
	cd ../Merced; $(MAKE)
	mv ../Merced/bin/merced .

upscatter: 
	cd ../fudge/processing/deterministic/upscatter; $(MAKE) CC=g++
	mv ../fudge/processing/deterministic/upscatter/bin/calcUpscatterKernel .

check:

clean:
	rm -f merced calcUpscatterKernel

realclean: clean
