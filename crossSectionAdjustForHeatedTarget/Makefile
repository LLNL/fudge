# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

PYTHON = python

dirList = Doc Src Test Python Misc include lib

.PHONY: lib all python check clean realclean cycleDirList

lib:
	cd Src; $(MAKE)
	cd include; $(MAKE) 
	cd lib; $(MAKE)

all:
	$(MAKE) dirList="Src Doc Python" cycleDirsTarget="all" cycleDirList

python:
	$(PYTHON) setup.py --quiet build

check:
	cd Test; $(MAKE) check

clean:
	rm -f *.pyc
	$(MAKE) cycleDirsTarget="clean" cycleDirList

realclean:
	rm -fr *.pyc */__pycache__ build
	$(MAKE) cycleDirsTarget="realclean" cycleDirList

cycleDirList:
	for d in $(dirList); do \
		cd $$d; \
		$(MAKE) $(cycleDirsTarget); \
		cd ..; \
	done
