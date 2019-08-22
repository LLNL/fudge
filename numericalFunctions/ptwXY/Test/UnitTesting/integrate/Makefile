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

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTWXY_PATH = ../../..
PTWXY_INCLUDE_PATH = $(PTWXY_PATH)/include
PTWXY_LIB_PATH = $(PTWXY_PATH)/lib
PTWXYUTILITIES_INCLUDE_PATH = $(PTWXY_PATH)/Test/Utilities
PTWXYUTILITIES_LIB_PATH = $(PTWXY_PATH)/Test/Utilities

CC = gcc
CDEBUGFLAGS = -g -O0 -Wall
ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTWXY_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -L$(PTWXY_LIB_PATH) \
    -I$(PTWXYUTILITIES_INCLUDE_PATH) -L$(PTWXYUTILITIES_LIB_PATH) -lptwXY_utilities -lnumericalfunctions
localCFLAGS = $(CDEBUGFLAGS) $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = groupTestAll1 integration integrationXY integrationWithWeight_x integrationWithWeight_sqrt_x integration2 runningIntegration groupTwoFunction

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

groupTestAll1: groupTestAll1.c
	$(CC) -o groupTestAll1 groupTestAll1.c -lptwXY -lm $(localCFLAGS)

integration: integration.c
	$(CC) -o integration integration.c -lptwXY -lm $(localCFLAGS)

integration2: integration2.c
	$(CC) -o integration2 integration2.c -lptwXY -lm $(localCFLAGS)

integrationXY: integrationXY.c
	$(CC) -o integrationXY integrationXY.c -lptwXY -lm $(localCFLAGS)

integrationWithWeight_x: integrationWithWeight_x.c
	$(CC) -o integrationWithWeight_x integrationWithWeight_x.c -lptwXY -lm $(localCFLAGS)

integrationWithWeight_sqrt_x: integrationWithWeight_sqrt_x.c
	$(CC) -o integrationWithWeight_sqrt_x integrationWithWeight_sqrt_x.c -lptwXY -lm $(localCFLAGS)

runningIntegration: runningIntegration.c
	$(CC) -o runningIntegration runningIntegration.c -lptwXY -lm $(localCFLAGS)

groupTwoFunction: groupTwoFunction.c
	$(CC) -o groupTwoFunction groupTwoFunction.c -lptwXY -lm $(localCFLAGS)

check: default
	for f in $(TARGETS); do $$f $(CHECKOPTIONS); done

clean:
	rm -f *.dat v

realclean: clean
	rm -f $(TARGETS)
