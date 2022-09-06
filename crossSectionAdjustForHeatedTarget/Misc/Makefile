# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

NUCLEAR_DIR = /usr/apps/nuclear/current/default
local_CFLAGS = $(CFLAGS) -I../Src  -I$(NUCLEAR_DIR)/include -g -O0

Objs = ../Src/crossSectionAdjustForHeatedTarget.o
LIBs = -L$(NUCLEAR_DIR)/lib -lnuclear

.PHONY: all clean realclean

all: nuc_xsec_adjust_for_heated_target_at_E nuc_xsec_adjust_for_heated_target_endl nuc_xsec_adjust_for_heated_target

nuc_xsec_adjust_for_heated_target_at_E: nuc_xsec_adjust_for_heated_target_at_E.c $(Objs)
	$(CC) $(local_CFLAGS) -o nuc_xsec_adjust_for_heated_target_at_E nuc_xsec_adjust_for_heated_target_at_E.c $(Objs) -lm

nuc_xsec_adjust_for_heated_target: nuc_xsec_adjust_for_heated_target.c $(Objs)
	$(CC) $(local_CFLAGS) -o nuc_xsec_adjust_for_heated_target nuc_xsec_adjust_for_heated_target.c $(Objs) -lm

nuc_xsec_adjust_for_heated_target_endl: nuc_xsec_adjust_for_heated_target_endl.c $(Objs)
	$(CC) $(local_CFLAGS) -o nuc_xsec_adjust_for_heated_target_endl nuc_xsec_adjust_for_heated_target_endl.c $(Objs) $(LIBs) -lm

clean:
	rm -f nuc_xsec_adjust_for_heated_target_endl nuc_xsec_adjust_for_heated_target_at_E nuc_xsec_adjust_for_heated_target

realclean: clean
