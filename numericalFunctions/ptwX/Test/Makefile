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

.PHONY: default check clean realclean doTarget

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTW1D_PATH = ..
PTW1D_INCLUDE_PATH = $(PTW1D_PATH)/include
PTW1D_LIB_PATH = $(PTW1D_PATH)/lib
PTW1D_ARCHIVE_FILE = $(PTW1D_LIB_PATH)/libptwX.a

CDEBUGFLAGS = -g -O0 -Wall
ptwXCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTW1D_INCLUDE_PATH) -L$(PTW1D_LIB_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lptwX -lnumericalfunctions -lm
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXCFLAGS)

subDirectories = 

default: testCreate

testCreate: testCreate.c $(PTWXY_ARCHIVE_FILE)
	$(CC) -o testCreate testCreate.c $(localCFLAGS)

check: default
	./testCreate $(CHECKOPTIONS)

clean:
	rm -f testCreate
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
