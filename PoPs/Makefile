# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SUB_DIRS = quantities groups families decays atomic Test

.PHONY: clean realclean doTarget

default:

clean:
	rm -f *.pyc
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	SAVED_PWD=`pwd`; \
	for directory in $(SUB_DIRS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done
