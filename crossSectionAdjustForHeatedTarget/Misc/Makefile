# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
