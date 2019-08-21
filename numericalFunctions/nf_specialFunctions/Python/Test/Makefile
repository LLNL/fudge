# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

subDirectories = UnitTesting

.PHONY: check checke clean realclean doTarget

check:
	$(MAKE) doTarget MODE=check

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done

