# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check clean realclean doTarget

subDirectories = UnitTesting

default:
	$(MAKE) doTarget MODE=''

check:
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) doTarget MODE=check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
