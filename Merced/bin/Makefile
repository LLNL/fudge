# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

.PHONY: default merced check clean realclean

default: merced

merced: 
	cd ../Src; $(MAKE)
	cp ../Src/merced .

check:

clean:
	rm -fr merced

realclean: clean
