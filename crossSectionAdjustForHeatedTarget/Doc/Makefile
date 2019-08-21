# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: all pdf crossSectionAdjustForHeatedTarget.pdf clean realclean

all: pdf

pdf: crossSectionAdjustForHeatedTarget.pdf

crossSectionAdjustForHeatedTarget.pdf: crossSectionAdjustForHeatedTarget.tex
	pdflatex crossSectionAdjustForHeatedTarget.tex

clean:
	rm -f *.log *.aux *.dvi *.pdf

realclean: clean
