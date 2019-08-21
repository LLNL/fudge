# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: doc dvi ps pdf clean realclean ptwXY.dvi

doc: pdf

dvi: ptwXY.dvi
ptwXY.dvi: ptwXY.tex
	latex ptwXY.tex
	latex ptwXY.tex
	latex ptwXY.tex

ps: ptwXY.ps
ptwXY.ps: ptwXY.dvi
	dvips -o ptwXY.ps ptwXY.dvi

pdf: ptwXY.pdf
ptwXY.pdf: ptwXY.dvi
	dvipdf ptwXY.dvi ptwXY.pdf

clean:
	rm -f *.dvi *.ps *.pdf *.toc *.lot
	cd Examples; $(MAKE) clean

realclean: clean
	rm -f *.aux *.log
	cd Examples; $(MAKE) realclean

