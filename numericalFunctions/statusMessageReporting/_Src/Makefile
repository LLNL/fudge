# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default collectSources clean realclean

default:

collectSources: clean
	cp ../Src/*.[ch] .

clean:
	rm -f *.c* *.h*

realclean: clean
