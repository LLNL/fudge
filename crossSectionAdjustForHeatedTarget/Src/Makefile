# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>


.PHONY: clean realclean

default: libcrossSectionAdjustForHeatedTarget.a

all: default

libcrossSectionAdjustForHeatedTarget.a: crossSectionAdjustForHeatedTarget.o
	$(AR) -cr libcrossSectionAdjustForHeatedTarget.a crossSectionAdjustForHeatedTarget.o

crossSectionAdjustForHeatedTarget.o: crossSectionAdjustForHeatedTarget.c
	$(CC) $(CFLAGS) -I. -c -o crossSectionAdjustForHeatedTarget.o crossSectionAdjustForHeatedTarget.c

clean:
	rm -f *.o *.a

realclean: clean
