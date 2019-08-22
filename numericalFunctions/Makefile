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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>

SOURCE_DIRS = ptwC ptwX ptwXY nf_Legendre nf_specialFunctions nf_integration
SUB_DIRS = include lib $(SOURCE_DIRS) Doc _Src

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

checkPython:
	$(MAKE) doTarget MODE=checkPython SUB_DIRS="$(SOURCE_DIRS)"

checkePython:
	$(MAKE) doTarget MODE=checkPython SUB_DIRS="$(SOURCE_DIRS)" CHECKOPTIONS=-e

collectSources:
	cd _Src; $(MAKE) collectSources

clean:
	rm -rf build
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf build *.pyc
	$(MAKE) doTarget MODE=realclean

doTarget:
	SAVED_PWD=`pwd`; \
	for directory in $(SUB_DIRS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done
