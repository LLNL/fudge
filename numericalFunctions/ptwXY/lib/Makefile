# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

sourceDir = ../Src

.PHONY: default clean realclean

default: libptwXY.a

libptwXY.a: $(sourceDir)/*.c $(sourceDir)/*.h
	cd $(sourceDir); $(MAKE)
	$(AR) -cr libptwXY.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
