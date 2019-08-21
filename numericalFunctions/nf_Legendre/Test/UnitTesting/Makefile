# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

subDirectories = GaussianQuadrature to_ptwXY normalize

.PHONY: default check clean realclean doTarget

default:
	$(MAKE) doTarget MODE=''

check: default
	-$(MAKE) doTarget MODE=check

clean:
	$(MAKE) doTarget MODE=clean

realclean: clean
	$(MAKE) doTarget MODE=realclean

doTarget:
	for directory in $(subDirectories); do cd $$directory; $(MAKE) $(MODE); cd ../; done
