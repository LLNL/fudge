# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PYTHON = python

subDirectories = initial UnitTesting Flat

.PHONY: check checke output clean realclean doTarget

TestFiles = t1.py  t2.py  t3.py  t4.py  t5.py  t6.py  t7.py  t8.py  t9.py

check:
	$(PYTHON) t7.py > /dev/null
	diff Temp/t7.*.out; if [ $$? != 0 ]; then echo "FAILURE: t7.py"; fi
	diff Temp/t7.orig.out Results/t7.orig.out; if [ $$? != 0 ]; then echo "FAILURE: t7.py"; fi
	$(PYTHON) t6.py > /dev/null
	diff Temp/t6.out Results/t6.out; if [ $$? != 0 ]; then echo "FAILURE: t6.py"; fi
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) check CHECKOPTIONS=-e

output:
	cd Temp; $(MAKE) realclean
	for f in $(TestFiles); do echo $$f; $$f > Temp/$$f.stdout 2> Temp/$$f.stderr; done
	cd Temp; $(PYTHON) output.py

clean:
	$(MAKE) doTarget MODE=clean
	cd Temp; $(MAKE) clean
	cd Utilities; $(MAKE) clean

realclean:
	$(MAKE) doTarget MODE=realclean
	cd Temp; $(MAKE) realclean
	cd Utilities; $(MAKE) realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
