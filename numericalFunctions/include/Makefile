# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default headers clean realclean

default: headers

headers:
	cp $(NFUPDATING_PATH)/include/*.h .

clean:
	rm -f *.h

realclean: clean
