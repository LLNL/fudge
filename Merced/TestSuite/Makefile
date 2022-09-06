# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

subDirs = Compton ENDFLegendre ENDFdoubleDiff GeneralEvaporation Kalbach Legendre Legendre2Body Madland Maxwell Watt captureGamma coherent doubleDifferential evaporation isotropic phaseSpace twoBody two_step uncorrelated weights
.PHONY: default clean realclean check

default:

clean:

realclean:
	for subdir in $(subDirs); do rm -f $$subdir/*.info 2> /dev/null; rm -f $$subdir/out*_new 2> /dev/null; done

check:
	./doTests.sh
