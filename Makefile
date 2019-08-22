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

.PHONY: clean realclean build extensions

EXTENSIONS = numericalFunctions crossSectionAdjustForHeatedTarget

default: inplace

build:
	python setup.py --quiet build

inplace:
	python setup.py --quiet build_ext --inplace build

clean:
	rm -rf build
	find . -name "*.pyc" -exec rm -f {} \;
	find . -name "*.so" -exec rm -f {} \;
	$(MAKE) extensions MODE=realclean

realclean: clean

extensions:
	SAVED_PWD=`pwd`; \
	for directory in $(EXTENSIONS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done

docs:
	epydoc-2.6 -o doc/code-ref --html --show-imports --exclude="fudge.processing.resonances.setup" --exclude="fudge.vis.gnuplot.fudge2dMultiPlot" --exclude="fudge.vis.gnuplot.endl[2-4]dplot" --no-frames --docformat="epytext" -v fudge
# --check --debug --graph="umlclasstree"

dist:
	python setup.py sdist --formats=gztar,zip

check: check-pqu check-nf check-fudge

check-heat: # BROKEN?
	python crossSectionAdjustForHeatedTarget/Python/Test/t.py

check-pqu:
	python pqu/test/testPhysicalQuantityWithUncertainty.py

check-nf:
	cd numericalFunctions/nf_specialFunctions/Python/Test/UnitTesting/; $(MAKE) check

FUDGETESTFILES = \
    fudge/core/math/test/testFudgeMath.py \
    fudge/core/math/test/test_linearAlgebra.py \
    fudge/processing/resonances/test/test_fudgeReconstructResonances.py  \
    fudge/processing/resonances/test/test_getScatteringMatrices.py \
    fudge/legacy/endl/test/test_endlProject.py \
    fudge/gnd/productData/distributions/test/__init__.py \
    fudge/gnd/reactionData/test/test_crossSection.py \
    fudge/gnd/test/testCovariances.py 
#    fudge/core/utilities/test/testUtilities.py \ # BROKEN?
#    fudge/core/math/test/testXYs.py \ # BROKEN?

check-fudge:
	for testFile in $(FUDGETESTFILES); do echo ; echo ======================================================================= ; echo \>\>\> TESTING $$testFile ; echo =======================================================================; echo ; python $$testFile; done
