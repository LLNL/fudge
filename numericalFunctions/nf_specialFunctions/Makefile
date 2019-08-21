# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default check checke clean realclean doTarget

subDirectories = Src Test lib include

default:
	cd lib; $(MAKE)
	cd include; $(MAKE)

check: default
	cd Test; $(MAKE) check

checke:
	$(MAKE) check CHECKOPTIONS=-e

clean:
	$(MAKE) doTarget MODE=clean

realclean:
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
