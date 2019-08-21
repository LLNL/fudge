# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

subDirectories = Utilities initial UnitTesting Flat

.PHONY: default check clean realclean doTarget

default:
	$(MAKE) doTarget MODE=''

check: default
	-$(MAKE) doTarget MODE=check

checke:
	$(MAKE) doTarget MODE=check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
