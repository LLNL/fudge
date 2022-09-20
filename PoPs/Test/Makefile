# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

PYTHON = python3

TARGETS = ids.py quantities.py gaugeBoson.py lepton.py baryon.py nuclearLevel.py isotope.py \
    chemicalElement.py database.py database2.py database3.py database4.py alias.py copyDBs.py \
	build.py iterate.py

default:

check:
	rm -rf Outputs
	mkdir Outputs
	for target in ${TARGETS}; do $(PYTHON) $$target > Outputs/$$target.out; ./check.com $$target $(CHECKOPTIONS); done

clean:
	rm -rf Outputs

realclean: clean
