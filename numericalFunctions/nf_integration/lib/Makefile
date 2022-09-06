# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

sourceDir = ../Src

.PHONY: clean realclean

libnf_integration.a: $(sourceDir)/*.c $(sourceDir)/*.h
	cd $(sourceDir); $(MAKE)
	$(AR) -cr libnf_integration.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
