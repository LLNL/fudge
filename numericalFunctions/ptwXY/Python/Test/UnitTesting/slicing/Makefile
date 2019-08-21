# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check clean realclean

default:

check:
	python getSlice.py
	python setSlice.py
	python xSlice.py

checke:
	$(MAKE) CHECKOPTIONS="-v" check

clean:

realclean: clean
