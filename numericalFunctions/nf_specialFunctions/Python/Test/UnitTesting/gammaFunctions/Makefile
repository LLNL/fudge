# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check clean realclean

default:

check:
	python gammaTest.py
	python incompleteGammaTest.py

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:

realclean: clean
