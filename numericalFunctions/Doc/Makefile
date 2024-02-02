# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: doc pdf clean realclean ptwXY.pdf

doc: pdf

pdf: ptwXY.pdf

ptwXY.pdf: ptwXY.tex
	pdflatex ptwXY.tex > /dev/null
	pdflatex ptwXY.tex > /dev/null
	pdflatex ptwXY.tex

clean:
	rm -f *.pdf *.toc *.lot
	cd Examples; $(MAKE) clean

realclean: clean
	rm -f *.aux *.log
	cd Examples; $(MAKE) realclean
