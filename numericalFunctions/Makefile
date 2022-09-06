# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python

SOURCE_DIRS = ptwC ptwX ptwXY nf_Legendre nf_specialFunctions nf_integration
SUB_DIRS = include lib $(SOURCE_DIRS) Doc _Src

.PHONY: default purgeIncludeAndLib includeAndLib $(SOURCE_DIRS) python utilities docs check checke checkPython checkePython \
	collectSources clean realclean doTarget statusMessageReporting

SMR_PATH = $(abspath ./statusMessageReporting)
SMR_INCLUDE_PATH = $(SMR_PATH)/include
SMR_LIB_PATH = $(SMR_PATH)/lib

MAKEFLAGS += SMR_PATH=$(SMR_PATH)

default: purgeIncludeAndLib $(SOURCE_DIRS)
	$(MAKE) utilities

purgeIncludeAndLib: 
	cd include; $(MAKE) realclean
	cd lib; $(MAKE) realclean

includeAndLib:
	cd include; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)
	cd lib; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)

statusMessageReporting:
	cd statusMessageReporting; $(MAKE)
	cp statusMessageReporting/include/*.h include
	cp statusMessageReporting/lib/*.a lib

ptwC: statusMessageReporting
	cd ptwC; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwC

ptwX: ptwC
	cd ptwX; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwX

ptwXY: ptwX
	cd ptwXY; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwXY

nf_Legendre: ptwXY
	cd nf_Legendre; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_Legendre

nf_specialFunctions: nf_Legendre
	cd nf_specialFunctions; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_specialFunctions

nf_integration: nf_specialFunctions                         # The chained dependencies are needed for the -j option. Why, do not know but does not work otherwise.
	cd nf_integration; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_integration

python:
	$(PYTHON) setup.py --quiet build

utilities:
	$(MAKE) doTarget MODE=utilities SUB_DIRS="$(SOURCE_DIRS)"

docs:

check:
	$(MAKE) doTarget MODE=check SUB_DIRS="$(SOURCE_DIRS)"

checke:
	$(MAKE) check CHECKOPTIONS=-e

checkPython:
	$(MAKE) doTarget MODE=checkPython SUB_DIRS="$(SOURCE_DIRS)"

checkePython:
	$(MAKE) CHECKOPTIONS=-e checkPython

collectSources:
	cd _Src; $(MAKE) collectSources

clean:
	rm -rf build
	cd statusMessageReporting; $(MAKE) clean
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf build *.pyc __pycache__
	cd statusMessageReporting; $(MAKE) realclean
	$(MAKE) doTarget MODE=realclean

doTarget: 
	SAVED_PWD=`pwd`; \
	for directory in $(SUB_DIRS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done
