# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>

.PHONY: default build inplace bin clean realclean build extensions

EXTENSIONS = statusMessageReporting numericalFunctions crossSectionAdjustForHeatedTarget Merced bin

EPYDOC=epydoc-2.6
PYTHON=python

default: inplace

build:
	$(PYTHON) setup.py --quiet build

inplace:
	$(PYTHON) setup.py --quiet build_ext --inplace build
	@echo 'INFO: This target does not build merced'

bin:
	cd bin; $(MAKE)

merced: bin

clean:
	rm -rf build
	find . -name "*.pyc" -exec rm -f {} \;
	find . -name "*.so" -exec rm -f {} \;
	$(MAKE) extensions MODE=realclean

realclean: clean
	rm -rf Out
	cd fudge; $(MAKE) realclean

pyclean:
	find . -name "*.pyc" -exec rm -f {} \;

extensions:
	SAVED_PWD=`pwd`; \
	for directory in $(EXTENSIONS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done

docs:
	cd doc/sphinx; $(MAKE) $(MODE) html; cd ../..

rebuild-test-data:
	cd fudge/gnds/covariances/test; python rebuild_test_data.py
	cd fudge/processing/resonances/test; python rebuild_test_data.py

dist:
	$(PYTHON) setup.py sdist --formats=gztar,zip

check: rebuild-test-data check-pqu check-nf check-smr check-fudge check-PoPs

check-site-packages:
	cd site_packages; $(MAKE) check

check-heat: # BROKEN?
	$(PYTHON) crossSectionAdjustForHeatedTarget/Python/Test/t.py

check-pqu:
	$(PYTHON) pqu/Check/check.py

check-PoPs:
	cd PoPs/Test/; $(MAKE) check

check-nf:
	cd numericalFunctions/nf_specialFunctions/Python/Test/UnitTesting/; $(MAKE) check

check-smr:
	cd statusMessageReporting/; $(MAKE)
	cd statusMessageReporting/Test/; $(MAKE) check

FUDGETESTFILES = \
    fudge/core/math/test/testFudgeMath.py \
    fudge/core/math/test/test_linearAlgebra.py \
    fudge/processing/resonances/test/test_reconstructResonances.py  \
    fudge/processing/resonances/test/test_getScatteringMatrices.py \
    fudge/legacy/endl/test/test_endlProject.py \
    fudge/gnds/productData/distributions/test/__init__.py \
    fudge/gnds/reactionData/test/test_crossSection.py \
    fudge/gnds/covariances/test/test_base.py \
    fudge/gnds/covariances/test/test_mixed.py \
    fudge/gnds/covariances/test/test_summed.py \
    fudge/gnds/covariances/test/test_covarianceSuite.py
#    fudge/particles/test/testParticles.py

FUDGECOVTESTFILES = \
    fudge/gnds/reactionData/test/test_crossSection.py \
    fudge/gnds/covariances/test/test_base.py \
    fudge/gnds/covariances/test/test_mixed.py \
    fudge/gnds/covariances/test/test_summed.py \
    fudge/gnds/covariances/test/test_covarianceSuite.py

check-cov:
	for testFile in $(FUDGECOVTESTFILES); do echo ; \
		echo ======================================================================= ; \
		echo \>\>\> TESTING $$testFile ; \
		echo =======================================================================; \
		echo ; python $$testFile; done

check-fudge:
	for testFile in $(FUDGETESTFILES); do echo ; \
		echo ======================================================================= ; \
		echo \>\>\> TESTING $$testFile ; \
		echo =======================================================================; \
		echo ; python $$testFile; done
