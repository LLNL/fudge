# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: doc dvi ps pdf clean realclean statusMessageReporting.dvi

doc: pdf

pdf: statusMessageReporting.pdf
statusMessageReporting.pdf: statusMessageReporting.dvi
	pdflatex statusMessageReporting.tex
	pdflatex statusMessageReporting.tex
	pdflatex statusMessageReporting.tex

dvi: statusMessageReporting.dvi
statusMessageReporting.dvi: statusMessageReporting.tex
	pdflatex statusMessageReporting.tex
	pdflatex statusMessageReporting.tex
	pdflatex statusMessageReporting.tex

ps: statusMessageReporting.ps
statusMessageReporting.ps:
	dvips -o statusMessageReporting.ps statusMessageReporting.dvi

clean:
	rm -f *.dvi *.ps *.pdf *.toc *.lot

realclean: clean
	rm -f *.aux *.log
