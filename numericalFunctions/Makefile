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

SOURCE_DIRS = ptwC ptwX ptwXY nf_Legendre nf_specialFunctions nf_integration
SUB_DIRS = include lib $(SOURCE_DIRS) Doc

.PHONY: default purgeIncludeAndLib includeAndLib $(SOURCE_DIRS) check checke clean realclean doTarget

SMR_PATH = $(abspath ../statusMessageReporting)
SMR_INCLUDE_PATH = $(SMR_PATH)/include
SMR_LIB_PATH = $(SMR_PATH)/lib

MAKEFLAGS += SMR_PATH=$(SMR_PATH)

default: purgeIncludeAndLib $(SOURCE_DIRS)

purgeIncludeAndLib:
	cd include; $(MAKE) realclean
	cd lib; $(MAKE) realclean

includeAndLib:
	cd include; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)
	cd lib; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)

ptwC:
	cd ptwC; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwC

ptwX: ptwC
	cd ptwX; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwX

ptwXY: ptwX
	cd ptwXY; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwXY

nf_Legendre: ptwX ptwXY
	cd nf_Legendre; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_Legendre

nf_specialFunctions: ptwC
	cd nf_specialFunctions; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_specialFunctions

nf_integration:
	cd nf_integration; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_integration

python:
	python setup.py --quiet build

check:
	$(MAKE) doTarget MODE=check SUB_DIRS="$(SOURCE_DIRS)"

checke:
	$(MAKE) doTarget MODE=check SUB_DIRS="$(SOURCE_DIRS)" CHECKOPTIONS=-e

clean:
	rm -rf build
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf build *.pyc
	$(MAKE) doTarget MODE=realclean

doTarget:
	SAVED_PWD=`pwd`; \
	for directory in $(SUB_DIRS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done
