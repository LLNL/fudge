# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

default: calcUpscatterKernel

calcUpscatterKernel: calcUpscatterKernel.c
	$(CC) $(CFLAGS) -g -o calcUpscatterKernel calcUpscatterKernel.c -lm

clean:
	rm -f calcUpscatterKernel
	rm -rf *.dSYM

realclean: clean
