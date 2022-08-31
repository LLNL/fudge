# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

coverageFolder = $(shell pwd)
srcDir = $(coverageFolder)/../Src
testDir = $(coverageFolder)/../TestSuite
coverageEXE = mercedWithCoverage
coveragePath = $(coverageFolder)/$(coverageEXE)
coverageSummary = coverage.info
subDirs = Compton ENDFLegendre ENDFdoubleDiff GeneralEvaporation Kalbach Legendre Legendre2Body \
          Madland Maxwell Watt captureGamma coherent doubleDifferential evaporation isotropic \
          phaseSpace twoBody two_step uncorrelated weights

.PHONY: default clean build coverage

default: coverage

coverage: $(coverageEXE)
	rm -rf individualCoverage.* html $(coverageSummary)
	for subdir in $(subDirs); do cd $(testDir)/$$subdir;\
	  for file in in.*; do \
	    rm -f $(srcDir)/*.gcda;\
	    $(coveragePath) $$file;\
	    lcov -c -d $(srcDir) -o "$(coverageFolder)/individualCoverage.$${subdir}_$${file}";\
	    rm -f utfil $(srcDir)/*.gcda;\
	  done;\
	done
	coverageList=(individualCoverage.*); coverageList=`echo $${coverageList[@]} | sed 's/ / -a /g'`;\
	lcov --add-tracefile $$coverageList -o $(coverageSummary)
	genhtml coverage.info -o html

$(coveragePath): build

build: clean
	cd $(srcDir); $(MAKE) VERSION=COVERAGE coverage; mv $(coverageEXE) $(coveragePath)

clean:
	@$(RM) -rf $(coveragePath) individualCoverage* html $(coverageSummary)
	cd $(srcDir); $(MAKE) clean
