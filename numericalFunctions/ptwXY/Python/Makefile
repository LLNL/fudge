# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default clean realclean doTarget

subDirectories = Src Test

default:
	cd Src; $(MAKE)

check: default
	cd Test; $(MAKE) check

checke: default
	cd Test; $(MAKE) check CHECKOPTIONS=-e

clean:
	rm -rf ptwXY ptwC ptwX
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf ptwXY ptwC ptwX
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
