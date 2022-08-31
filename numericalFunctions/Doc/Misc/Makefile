# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

NUMERICALFUNCTIONS_PATH = ../../
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = infill_f3_f4

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

infill_f3_f4: infill_f3_f4.c
	$(CC) -o infill_f3_f4 infill_f3_f4.c $(localCFLAGS) -lm

check:

clean:
	rm -f $(TARGETS)

realclean: clean
