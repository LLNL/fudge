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

subDirectories = initial UnitTesting Flat

.PHONY: check checke output clean realclean doTarget

TestFiles = t1.py  t2.py  t3.py  t4.py  t5.py  t6.py  t7.py  t8.py  t9.py

check:
	./t7.py > /dev/null
	diff Temp/t7.*.out; if [ $$? != 0 ]; then echo "FAILURE: t7.py"; fi
	diff Temp/t7.orig.out Results/t7.orig.out; if [ $$? != 0 ]; then echo "FAILURE: t7.py"; fi
	./t6.py > /dev/null
	diff Temp/t6.out Results/t6.out; if [ $$? != 0 ]; then echo "FAILURE: t6.py"; fi
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) check CHECKOPTIONS=-e

output:
	cd Temp; $(MAKE) realclean
	for f in $(TestFiles); do echo $$f; $$f > Temp/$$f.stdout 2> Temp/$$f.stderr; done
	cd Temp; python output.py

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
